# Run the pooled-testing experiments, writing every solve to the SQLite store
# (data/solves.db). Analysis lives separately in analyse.jl and reads the store,
# so results can be analysed while a run is still in progress.
#
#     julia --project=. -t 4 run.jl                   # all experiments
#     julia --project=. -t 4 run.jl --experiments 1,3 # a subset
#
# Solves already in the store are skipped, so an interrupted/reclaimed run simply
# resumes by re-running the same command.

using PooledTesting, CSV, DataFrames, Distributions, Random, ArgParse

# Experiment parameters. Point POOLED_CONSTANTS at another file (e.g.
# constants_test.jl) to run the whole pipeline quickly with tiny parameters.
include(abspath(get(ENV, "POOLED_CONSTANTS", joinpath(@__DIR__, "constants.jl"))))

const ROOT = @__DIR__   # pilotdata.csv and the data/ store live alongside this script

"""
Approximate MILP warm-started from a fresh greedy solution, which also scales
the relative error budget `accuracy` when given (see `accuracy_params`).
Greedy takes ~1s against minutes for the MILP, so recomputing it here is
simpler than plumbing pools from the separate `:greedy` experiment solves.
"""
function warm_approximate(pop; T, G, verbose=false, kwargs...)
    _, start, _ = greedy(pop; T=T, G=G, verbose=verbose)
    return approximate(pop; T=T, G=G, start=start, verbose=verbose, kwargs...)
end

# The standard approx-vs-greedy algorithm pair. With an accuracy budget, the
# MILP derives (K, MIPGapAbs) per instance via accuracy_params; with
# accuracy=nothing it uses the fixed MILP_K and the global MIPGap.
approx_vs_greedy(accuracy=nothing) = [
    (name=:approx, fn=warm_approximate,
     args=accuracy === nothing ? Dict(:K => MILP_K) : Dict(:accuracy => accuracy)),
    (name=:greedy, fn=greedy, args=Dict()),
]

"Pilot population from the CSV, utilities scaled to integers in [1, UTIL_UPPER_BOUND]."
function pilot_population()
    df = DataFrame(CSV.File(joinpath(ROOT, "pilotdata.csv")))
    factor = UTIL_UPPER_BOUND / maximum(df.baseline_utility)
    pop = Population{Int}()
    for row in eachrow(df)
        pop[row.id] = (row.health_probability,
                       clamp(round(Int, row.baseline_utility*factor), 1, UTIL_UPPER_BOUND))
    end
    return pop
end

"`SYNTHETIC_REPS` populations of size `SYNTHETIC_N`, health ~ U(0.5,1), utilities ~ fitted Normal."
function synthetic_populations()
    df = DataFrame(CSV.File(joinpath(ROOT, "pilotdata.csv")))
    udist = fit(Normal, (UTIL_UPPER_BOUND/maximum(df.baseline_utility)) .* df.baseline_utility)
    return [generate_instance(SYNTHETIC_N, Uniform(0.5, 1), udist) for _ in 1:SYNTHETIC_REPS]
end

"Run experiment `num` (per EXPERIMENT_SPECS) into the store `db`."
function run_experiment(db, num)
    spec = EXPERIMENT_SPECS[num]
    Random.seed!(spec.seed)
    if spec.kind == :pilot
        run_experiments(db, approx_vs_greedy(spec.accuracy), [pilot_population()], PILOT_BUDGETS, [spec.G];
                        experiment=num, multithread=true)
    elseif spec.kind == :synthetic
        run_experiments(db, approx_vs_greedy(spec.accuracy), synthetic_populations(), SYNTHETIC_BUDGETS, [spec.G];
                        experiment=num, multithread=true)
    elseif spec.kind == :overlap
        algs = [(name=:disjoint, fn=exact, args=Dict(:k => 1)),
                (name=:two_overlap, fn=exact, args=Dict(:k => 2))]
        pops = [generate_instance(10, 0:0.1:1, 1:3) for _ in 1:SYNTHETIC_REPS]
        old = (mipgap=PooledTesting.GUROBI_MIPGAP[],
               threads=PooledTesting.GUROBI_THREADS[],
               params=copy(PooledTesting.GUROBI_EXTRA_PARAMS[]),
               param_file=PooledTesting.GUROBI_PARAM_FILE[])
        configure_gurobi!(mipgap=get(spec, :mipgap, old.mipgap), threads=old.threads,
                          params=old.params, param_file=old.param_file)
        try
            run_experiments(db, algs, pops, [2, 3, 4, 5], [spec.G]; experiment=num, multithread=true)
        finally
            configure_gurobi!(; old...)
        end
    end
end

"Warm up Julia so experiment timings aren't polluted by first-call compilation."
function warmup()
    println("\nWARMING UP THE ENGINE")
    pop = generate_instance(10, 0:0.1:1, 1:10)
    warm_approximate(pop; T=2, G=5, K=MILP_K)  # covers greedy + warm-started MILP
    println("READY TO RUMBLE\n")
end


## COMMAND LINE INTERFACE

const EXPERIMENT_NUMBERS = sort(collect(keys(EXPERIMENT_SPECS)))

function parse_commandline(args)
    s = ArgParseSettings(description="Run pooled-testing experiments into the solve store.")
    @add_arg_table! s begin
        "--experiments"
            help = "Comma-separated experiments to run. Available: $(join(EXPERIMENT_NUMBERS, ", ")). Default: all."
            arg_type = String
            default = join(EXPERIMENT_NUMBERS, ",")
        "--rootdir"
            help = "Directory for the solve store (default: this experiments directory)."
            arg_type = String
            default = ROOT
    end
    return parse_args(args, s)
end

function main(args)
    parsed = parse_commandline(args)
    db = open_store(parsed["rootdir"])
    warmup()
    for s in split(parsed["experiments"], [',', ' ']; keepempty=false)
        num = tryparse(Int, strip(s))
        if num === nothing || !haskey(EXPERIMENT_SPECS, num)
            @warn "Unknown experiment: $(s). Available: $(join(EXPERIMENT_NUMBERS, ", "))"
            continue
        end
        println("\nSTARTING EXPERIMENT $(num)")
        run_experiment(db, num)
        @info "Completed experiment $(num)"
    end
    @info "All requested experiments done. Analyse with: julia --project=. analyse.jl"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
