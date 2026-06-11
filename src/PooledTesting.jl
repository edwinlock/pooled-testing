"""
    PooledTesting

Welfare-maximising pooled testing: exact and approximate algorithms for computing
test allocations, and a SQLite-backed store of solve results.

The experiment scripts in `experiments/` (`run.jl`, `analyse.jl`) `using` this
package; see those for how the experiments are driven.
"""
module PooledTesting

using JuMP, Gurobi, MosekTools, MathOptInterface
using Roots, Combinatorics, DataStructures
using SQLite, DBInterface, DataFrames, SHA, Dates, Pkg
using ProgressMeter, PrecompileTools

# ---- Gurobi solver settings -------------------------------------------------
# Applied to every MILP via `set_gurobi_params!`. Defaults match the paper; tune
# with `configure_gurobi!`. MIPGap must stay tighter than the MILP-vs-greedy
# welfare difference or the comparison is invalid; threads are capped low because
# MILP branch-and-bound parallelises poorly past a few.
const GUROBI_MIPGAP  = Ref(1e-4)
const GUROBI_THREADS = Ref(8)

"Set the Gurobi MIPGap and/or thread count used by every MILP solve."
function configure_gurobi!(; mipgap=GUROBI_MIPGAP[], threads=GUROBI_THREADS[])
    GUROBI_MIPGAP[] = mipgap
    GUROBI_THREADS[] = threads
    return nothing
end

set_gurobi_params!(m) = set_optimizer_attributes(m,
    "Threads" => GUROBI_THREADS[], "MIPGap" => GUROBI_MIPGAP[])

# One Gurobi environment per thread (not safe to share across concurrent solves),
# created lazily. A Dict avoids assuming the thread count at precompile time.
const GRB_ENVS = Dict{Int, Gurobi.Env}()
const GRB_ENV_LOCK = ReentrantLock()
grb_env() = lock(GRB_ENV_LOCK) do
    get!(() -> Gurobi.Env(), GRB_ENVS, Threads.threadid())
end

# Gurobi environments are C handles that cannot survive precompilation, so discard
# any created during the precompile workload and start each session fresh.
function __init__()
    empty!(GRB_ENVS)
end

# A population maps participant ids to (health probability, integer utility).
const Population{T} = Dict{T, Tuple{Float64, Int}} where T

# ---- Library ----------------------------------------------------------------
include("models/utils.jl")
include("models/approximation-models.jl")
include("models/linear-models.jl")
include("models/conic-models.jl")
include("algorithms.jl")
include("datastore.jl")
include("runner.jl")

# ---- Public API -------------------------------------------------------------
export Population, configure_gurobi!
export greedy, exact, approximate, conic           # algorithms
export welfare, generate_instance, scale_utilities # population utilities
export milp_guarantee                              # approximation guarantee
export open_store, record_population!, record_solve!, solved_keys, load_solves,
       population_hash, param_key, store_path       # solve store
export run_experiments, store_solve!, main_thread_only  # runner

# ---- Precompilation ---------------------------------------------------------
# Bake the solve + store pipeline into the precompile cache so the first call in
# a session isn't paid as JIT latency. The solver calls are guarded: if no Gurobi
# /MOSEK licence is present, precompilation still succeeds (just covers less).
@setup_workload begin
    pop = Population{Int}(1 => (0.9, 5), 2 => (0.8, 3), 3 => (0.7, 4))
    @compile_workload begin
        milp_guarantee([0.9, 0.8], [5, 3]; T=2, G=2, K=5)
        db = open_store(mktempdir())
        h = record_population!(db, pop; experiment=1, pop_index=1)
        record_solve!(db, h, 2, 2, (name=:approx, fn=approximate, args=Dict(:K=>5)), 1.0, 0.1, 10)
        solved_keys(db); load_solves(db)
        try
            approximate(pop; T=1, G=2, K=5)   # Gurobi
            greedy(pop; T=1, G=2)             # MOSEK
        catch
        end
    end
end

end # module
