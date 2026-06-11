# Run the pooled-testing experiments, writing every solve to the SQLite store
# (data/solves.db). Analysis (tables, plots) lives separately in analyse.jl and
# reads the store — so results can be analysed while a run is still in progress.
#
#     julia --project=. -t 4 run.jl                      # all experiments
#     julia --project=. -t 4 run.jl --experiments 1,3    # a subset
#
# Solves already in the store are skipped, so an interrupted/reclaimed run simply
# resumes by re-running the same command.

using CSV, DataFrames, Statistics, Dates, Distributions, ProgressMeter,
      DataStructures, Gurobi, ArgParse, Random

# Gurobi solver settings (used by the MILP models in optimisation.jl).
const GUROBI_MIPGAP = 1e-4    # 0.01% — matches the paper's default; must stay tighter
                              # than the MILP-vs-greedy difference or the comparison is invalid.
const GUROBI_THREADS = 8      # MILP parallelises poorly past a few threads (benchmarked).

include("optimisation.jl")
include("datastore.jl")


## RUNNER

# MOSEK (via greedy's conic subroutine) aborts when called from a Julia worker
# thread on Linux aarch64, so greedy solves run on the main thread; Gurobi-based
# algorithms run in parallel.
main_thread_only(alg) = alg.fn === greedy

"Solve one work item `(pop_hash, pop, T, G, alg)` and record it in the store."
function store_solve!(db, (h, pop, T, G, alg))
    start = Dates.now()
    welfare, _, guarantee = alg.fn(pop; T=T, G=G, alg.args...)
    record_solve!(db, h, T, G, alg, welfare, guarantee, Dates.value(Dates.now() - start))
end

"""
Solve every (population, budget, poolsize, algorithm) combination not already in
the store, recording each result. Gurobi solves run in parallel (`multithread`);
MOSEK (greedy) solves always run serially on the main thread.
"""
function run_experiments(db, algs, populations, budgets, poolsizes; experiment, multithread=false)
    hashes = [record_population!(db, pop; experiment, pop_index=i)
              for (i, pop) in enumerate(populations)]
    done = solved_keys(db)
    todo(h, T, G, alg) = (h, T, G, String(alg.name), param_key(alg)) ∉ done

    work = [(h, pop, T, G, alg)
            for (h, pop) in zip(hashes, populations)
            for T in budgets for G in poolsizes for alg in algs if todo(h, T, G, alg)]
    nskip = length(hashes)*length(budgets)*length(poolsizes)*length(algs) - length(work)
    nskip > 0 && println("  ($(nskip) solves already in store, skipping)")

    progress = Progress(length(work); desc="Solving ", showspeed=true)
    desc(w) = "$(w[5].name) G=$(w[4]) B=$(w[3]) "
    par, serial = filter(w -> !main_thread_only(w[5]), work), filter(w -> main_thread_only(w[5]), work)

    if multithread && !isempty(par)
        # :greedy scheduling pulls cells one at a time, so no thread idles on a slow solve.
        Threads.@threads :greedy for w in par
            store_solve!(db, w); ProgressMeter.next!(progress; desc=desc(w))
        end
    else
        for w in par
            ProgressMeter.update!(progress; desc=desc(w)); store_solve!(db, w)
            ProgressMeter.next!(progress; desc=desc(w))
        end
    end
    for w in serial  # MOSEK: main thread only
        ProgressMeter.update!(progress; desc=desc(w)); store_solve!(db, w)
        ProgressMeter.next!(progress; desc=desc(w))
    end
    return nothing
end


## POPULATIONS AND ALGORITHMS

const UTIL_UPPER_BOUND = 50
const PILOT_BUDGETS = collect(2:4:30)      # {2, 6, ..., 30}
const SYNTHETIC_BUDGETS = collect(2:2:12)  # {2, 4, ..., 12}

"Extract the pilot population from a CSV, scaling utilities to integers in [1, upper]."
function extract_population(filename, upper)
    df = DataFrame(CSV.File(filename))
    factor = upper / maximum(df.baseline_utility)
    pop = Population{Int}()
    for row in eachrow(df)
        pop[row.id] = (row.health_probability, clamp(round(Int, row.baseline_utility*factor), 1, upper))
    end
    return pop
end

"""
`reps` synthetic populations of size `n`, with health probabilities from
`Uniform(0.5, 1)` and utilities from a Normal fitted to the pilot's utilities.
"""
function synthetic_populations(rootdir, n, reps)
    df = DataFrame(CSV.File(joinpath(rootdir, "pilotdata.csv")))
    udist = fit(Normal, (UTIL_UPPER_BOUND/maximum(df.baseline_utility)) .* df.baseline_utility)
    return [generate_instance(n, Uniform(0.5, 1), udist) for _ in 1:reps]
end

# Standard approx-vs-greedy algorithm pair, parameterised by approximation parameter K.
approx_vs_greedy(K) = [
    (name=:approx, fn=approximate, args=Dict(:K => K, :verbose => false)),
    (name=:greedy, fn=greedy, args=Dict()),
]


## EXPERIMENTS

"Approx-vs-greedy on the pilot population (experiments 1–2)."
function pilot_experiment(db, num; rootdir, G, K, seed)
    Random.seed!(seed)
    pop = extract_population(joinpath(rootdir, "pilotdata.csv"), UTIL_UPPER_BOUND)
    run_experiments(db, approx_vs_greedy(K), [pop], PILOT_BUDGETS, [G]; experiment=num, multithread=true)
end

"Approx-vs-greedy on `reps` synthetic populations of size `n` (experiments 3–4)."
function synthetic_experiment(db, num; rootdir, n, G, K, reps, seed)
    Random.seed!(seed)
    pops = synthetic_populations(rootdir, n, reps)
    run_experiments(db, approx_vs_greedy(K), pops, SYNTHETIC_BUDGETS, [G]; experiment=num, multithread=true)
end

"Disjoint vs 2-overlap on small populations (experiment 5)."
function overlap_experiment(db, num; rootdir, seed)
    Random.seed!(seed)
    pops = [generate_instance(10, 0:0.1:1, 1:3) for _ in 1:20]
    algs = [(name=:disjoint, fn=exact, args=Dict(:k => 1)),
            (name=:two_overlap, fn=exact, args=Dict(:k => 2))]
    run_experiments(db, algs, pops, [2, 3, 4, 5], [10]; experiment=num, multithread=true)
end

const EXPERIMENTS = Dict(
    "1" => (db, rootdir) -> pilot_experiment(db, 1; rootdir, G=5, K=15, seed=1002),
    "2" => (db, rootdir) -> pilot_experiment(db, 2; rootdir, G=10, K=15, seed=1003),
    "3" => (db, rootdir) -> synthetic_experiment(db, 3; rootdir, n=150, G=5, K=15, reps=20, seed=2001),
    "4" => (db, rootdir) -> synthetic_experiment(db, 4; rootdir, n=150, G=10, K=15, reps=20, seed=2002),
    "5" => (db, rootdir) -> overlap_experiment(db, 5; rootdir, seed=3001),
)
const EXPERIMENT_NUMBERS = sort(collect(keys(EXPERIMENTS)))

"Warm up Julia so experiment timings aren't polluted by first-call compilation."
function warmup()
    println("\nWARMING UP THE ENGINE")
    pop = generate_instance(10, 0:0.1:1, 1:10)
    greedy(pop; T=2, G=5); approximate(pop; T=2, G=5, K=15)
    println("READY TO RUMBLE\n")
end


## COMMAND LINE INTERFACE

function parse_commandline(args)
    s = ArgParseSettings(description="Run pooled-testing experiments into the solve store.")
    @add_arg_table! s begin
        "--experiments"
            help = "Comma-separated experiments to run. Available: $(join(EXPERIMENT_NUMBERS, ", ")). Default: all."
            arg_type = String
            default = join(EXPERIMENT_NUMBERS, ",")
        "--rootdir"
            help = "Root directory for input data and the solve store (default: current directory)."
            arg_type = String
            default = "."
    end
    return parse_args(args, s)
end

function main(args)
    parsed = parse_commandline(args)
    rootdir = parsed["rootdir"]
    requested = split(parsed["experiments"], [',', ' ']; keepempty=false)

    db = open_store(rootdir)
    warmup()

    for key in strip.(requested)
        if haskey(EXPERIMENTS, key)
            println("\nSTARTING EXPERIMENT $(key)")
            try
                EXPERIMENTS[key](db, rootdir)
                @info "Completed experiment $(key)"
            catch e
                @error "Failed to run experiment $(key)" exception=(e, catch_backtrace())
            end
        else
            @warn "Unknown experiment: $(key). Available: $(join(EXPERIMENT_NUMBERS, ", "))"
        end
    end
    @info "All requested experiments done. Analyse with: julia --project=. analyse.jl"
end

# Run only when executed as a script (not when included).
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
