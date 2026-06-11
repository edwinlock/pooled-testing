# Per-solve result store, backed by SQLite.
#
# Two tables:
#   * `populations` — one row per distinct population, keyed by a content hash of
#     its (q, u) data. Records its size, the experiment and seed it came from, and
#     the data itself, so any population is fully recoverable.
#   * `solves` — one row per (population, budget, poolsize, algorithm, params)
#     solve. Keyed independently of the experiment, so a solve is *reused* whenever
#     the same population/cell/params recur — across reruns and across experiments.
#
# This gives resumability finer than the per-experiment guard (an interrupted run
# loses at most the in-flight solves) and live analysis (`analyse.jl` can query
# the database while a run writes — SQLite WAL mode allows concurrent readers).
# SQLite serialises writes and commits atomically, so the store stays consistent
# under the parallel solve loop and survives a process kill mid-write.

# (SQLite, DBInterface, DataFrames, SHA, Dates, Pkg are `using`-ed by the module.)

# Provenance captured once per process. Set EC2_INSTANCE_TYPE in the environment
# to record which machine produced a solve.
const INSTANCE_TYPE = get(ENV, "EC2_INSTANCE_TYPE", "$(Sys.CPU_THREADS)core-$(Sys.ARCH)")
gurobi_version() = try
    string(Pkg.dependencies()[Base.UUID("2e9cd046-0924-5485-92f1-d5272153d98b")].version)
catch
    "unknown"
end
const GUROBI_VERSION = gurobi_version()

# Algorithm args that change cosmetic output but not the solution (excluded from
# the solve key, so e.g. a verbose and non-verbose solve are treated as identical).
const NON_SOLVE_ARGS = (:verbose,)

store_path(rootdir) = joinpath(rootdir, "data", "solves.db")

"""
Content hash of a population: a stable key independent of run order or how the
population was generated. Hashes the (q, u) pairs in canonical (sorted) order.
"""
population_hash(pop) = bytes2hex(sha1(string(sort([(float(v[1]), Int(v[2])) for v in values(pop)]))))[1:16]

"""
Canonical string of the solve-affecting parameters of an algorithm: its own args
(e.g. K for the MILP, k for overlap models, minus cosmetic ones) plus the global
MIPGap. Two solves with the same params are interchangeable.
"""
param_key(alg) = join(["$p=$v" for (p, v) in sort(collect(alg.args)) if p ∉ NON_SOLVE_ARGS] ∪
                      ["mipgap=$(GUROBI_MIPGAP[])"], ";")

"""
Open (creating if needed) the store and ensure its schema. WAL mode lets a reader
(analyse.jl) query while a run writes. Primary keys make both tables idempotent.
"""
function open_store(rootdir)
    mkpath(joinpath(rootdir, "data"))
    db = SQLite.DB(store_path(rootdir))
    DBInterface.execute(db, "PRAGMA journal_mode=WAL;")     # concurrent readers while writing
    DBInterface.execute(db, "PRAGMA busy_timeout=5000;")    # wait, don't error, on a brief lock
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS populations (
            pop_hash TEXT PRIMARY KEY,
            experiment INTEGER, pop_index INTEGER, n INTEGER,
            q TEXT, u TEXT
        );""")
    DBInterface.execute(db, """
        CREATE TABLE IF NOT EXISTS solves (
            pop_hash TEXT, budget INTEGER, poolsize INTEGER, alg TEXT, params TEXT,
            welfare REAL, guarantee REAL, time_ms INTEGER,
            instance_type TEXT, gurobi_version TEXT, timestamp TEXT,
            PRIMARY KEY (pop_hash, budget, poolsize, alg, params)
        );""")
    return db
end

"Record a population's metadata (idempotent on its hash). Returns the hash."
function record_population!(db, pop; experiment, pop_index)
    h = population_hash(pop)
    q = join((v[1] for v in values(pop)), ",")
    u = join((v[2] for v in values(pop)), ",")
    DBInterface.execute(db, """
        INSERT OR IGNORE INTO populations (pop_hash, experiment, pop_index, n, q, u)
        VALUES (?,?,?,?,?,?);""", (h, experiment, pop_index, length(pop), q, u))
    return h
end

"""
Set of solve keys `(pop_hash, budget, poolsize, alg, params)` already stored, so a
run can skip work already done.
"""
function solved_keys(db)
    rows = DBInterface.execute(db, "SELECT pop_hash, budget, poolsize, alg, params FROM solves") |> DataFrame
    return Set((r.pop_hash, r.budget, r.poolsize, r.alg, r.params) for r in eachrow(rows))
end

"Record one completed solve (atomic, idempotent on its key)."
function record_solve!(db, pop_hash, T, G, alg, welfare, guarantee, time_ms)
    DBInterface.execute(db, """
        INSERT OR REPLACE INTO solves
        (pop_hash, budget, poolsize, alg, params, welfare, guarantee, time_ms,
         instance_type, gurobi_version, timestamp)
        VALUES (?,?,?,?,?,?,?,?,?,?,?);""",
        (pop_hash, T, G, String(alg.name), param_key(alg), welfare, guarantee,
         Int(round(time_ms)), INSTANCE_TYPE, GUROBI_VERSION, string(Dates.now())))
    return nothing
end

"Load `solves` joined with `populations` as a DataFrame, for analysis."
load_solves(db::SQLite.DB) = DBInterface.execute(db,
    "SELECT s.*, p.experiment, p.pop_index, p.n
     FROM solves s LEFT JOIN populations p ON s.pop_hash = p.pop_hash") |> DataFrame
load_solves(rootdir::AbstractString) = load_solves(open_store(rootdir))
