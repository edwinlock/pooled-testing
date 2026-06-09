# Pooled testing experiments.
#
# Run from the command line, selecting which experiments to run:
#
#     julia experiments.jl --experiments 1,3,5
#     julia experiments.jl                       # runs all experiments
#     julia experiments.jl --rootdir path/to/dir # set output root directory
#
# Each experiment is a self-contained `experimentN(; rootdir)` function.

# Load packages
using CSV, DataFrames, Statistics, Dates, StatsPlots, Distributions, ProgressMeter, DataStructures, Gurobi, ArgParse, Random

# Include optimisation code
include("optimisation.jl")


## HELPERS

"""
Solve a single (population, budget, poolsize) cell with every algorithm and
return the result row as a `Dict`. Each call is independent of every other, so
this is the unit of work parallelised by `run_experiments`.
"""
function run_cell(algs, i, pop, T, G)
    result = Dict{Symbol, Any}(:budget => T, :population => i, :poolsize => G)
    for (name, fn, args) in algs
        start = Dates.now()
        w, pools, error = fn(pop; T=T, G=G, args...)
        time = Dates.now() - start
        result[Symbol("$(name)_welfare")] = w
        result[Symbol("$(name)_error")] = error
        result[Symbol("$(name)_time")] = time
    end
    return result
end


"""
Convenience function to run experiments with specified algorithms, populations, budgets and poolsizes.

Set `multithread=true` to solve the independent (population, budget, poolsize)
cells in parallel across Julia threads (start Julia with e.g. `julia -t auto`).
Recorded solve times stay per-solve wall-clock, but for clean timings cap the
solver's own threads (e.g. Gurobi `Threads`=1) so concurrent solves don't
contend, and keep the timed experiments single-threaded.
"""
function run_experiments(algs, populations, budgets, poolsizes; multithread=false)
    df = DataFrame(:budget=>Int[], :population=>Int[], :poolsize => Int[])  # for storing results
    for alg in algs  # add columns for all algorithms
        df[!, "$(alg.name)_welfare"] = Float64[]
        df[!, "$(alg.name)_error"] = Float64[]
        df[!, "$(alg.name)_time"] = Millisecond[]
    end
    # Flatten the loops into a flat list of independent work items.
    work = [(i, pop, T, G)
            for (i, pop) in enumerate(populations) for T in budgets for G in poolsizes]
    results = Vector{Dict{Symbol, Any}}(undef, length(work))
    progress = Progress(length(work); desc="Running experiments ")  # ProgressMeter is thread-safe
    # Current cell shown inline in the bar's description (budget, pool size).
    celldesc(T, G) = "Running G=$(G), B=$(T) "
    if multithread
        # Each thread writes its own index, so no synchronisation is needed.
        Threads.@threads for k in eachindex(work)
            i, pop, T, G = work[k]
            results[k] = run_cell(algs, i, pop, T, G)
            ProgressMeter.next!(progress; desc=celldesc(T, G))
        end
    else
        for k in eachindex(work)
            i, pop, T, G = work[k]
            # Update the description to the cell about to run, so a long solve
            # shows its own budget/G in the bar while it churns.
            ProgressMeter.update!(progress; desc=celldesc(T, G))
            results[k] = run_cell(algs, i, pop, T, G)
            ProgressMeter.next!(progress; desc=celldesc(T, G))
        end
    end
    # Assemble the DataFrame serially (DataFrame push! is not thread-safe).
    for result in results
        push!(df, result)
    end
    return df
end


"""
Compute comparisons between the two algorithms specified in `algs`.
"""
function add_comparisons!(df, algs)
    symb1, symb2 = ["$(alg.name)_welfare" for alg in algs]
    # Add diffs and ratios to df
    df[!, :diff] = df[!,symb1] - df[!,symb2]
    df[!, :ratio] = df[!,symb1] ./ df[!,symb2]
end


"""Compute summary table from df."""
# Mean of a collection of `Millisecond` durations, rounded to whole milliseconds.
# `mean` on `Millisecond` returns a fractional period, so round and rebuild a
# `Millisecond` to keep the column a clean duration type.
roundmean(x) = Millisecond(round(Int, Dates.value(sum(x)) / length(x)))
function create_summary(df)
    output = combine(
        groupby(df, :budget),
        :approx_welfare => mean,
        :approx_error => mean => "Add. error",
        :approx_time => roundmean => "approx time",
        :greedy_welfare => mean,
        :greedy_time => roundmean => "greedy time",
    )
    return output
end


"""Plot welfares of approx and greedy algorithm from df."""
function plot_welfares(df, filename)
    @df df violin(
        :budget,
        :approx_welfare,
        side=:left,
        linewidth=0,
        label="",
        legend=false,
        xticks=2:2:12,
        xlabel="Test budget", ylabel="Welfare",
        ylims=nothing,
        size=(400,300))
    @df df dotplot!(:budget, :approx_welfare, side=:left, marker=(:black,stroke(0)), label="")
    @df df violin!(:budget, :greedy_welfare, side=:right, linewidth=0, label="Greedy")
    @df df dotplot!(:budget, :greedy_welfare, side=:right, marker=(:black,stroke(0)), label="")
    Plots.pdf(filename)
end


"""Plot welfare ratios of approx and greedy algorithm from df."""
function plot_ratios(df, filename)
    @df df violin(
        :budget,
        :ratio,
        linewidth=0,
        label="",
        legend=:right,
        xlabel="Test budget",
        ylabel="Welfare ratio",
        xticks=2:2:12,
        ylims=nothing,
        size=(400,300)
        )
    @df df dotplot!(:budget, :ratio, marker=(:black,stroke(0)), label="")
    Plots.pdf(filename)
end


"""
Scatter plot of disjoint welfare (x) against two-overlap welfare (y), with a
different colour and marker shape for each test budget, plus a y=x reference
line. An alternative to the violin plot for comparing the two allocations.
"""
function plot_disjoint_vs_overlap(df, filename)
    @df df scatter(:disjoint_welfare, :two_overlap_welfare,
        group=:budget,
        markershape=:auto,    # distinct shape per budget group
        xlabel="Disjoint welfare",
        ylabel="2-overlap welfare",
        legend=:topleft,
        legendtitle="Test budget",
        alpha=0.7,
        size=(400,300))
    # Diagonal y=x reference line: points above it are where 2-overlap wins.
    lo = min(minimum(df.disjoint_welfare), minimum(df.two_overlap_welfare))
    hi = max(maximum(df.disjoint_welfare), maximum(df.two_overlap_welfare))
    plot!([lo, hi], [lo, hi], line=:dash, color=:gray, label="y=x")
    Plots.pdf(filename)
end


"Extract population data from CSV file"
function extract_population(filename, utility_upper_bound)
    df = DataFrame(CSV.File(filename))
    scaling_factor = utility_upper_bound / maximum(df.baseline_utility)
    trial_population = Population{Int}()  # with integral utilities
    for row in eachrow(df)
        q = row.health_probability
        u = clamp(round(Int, row.baseline_utility*scaling_factor), 1, utility_upper_bound)
        trial_population[row.id] = (q, u)
    end
    return trial_population, scaling_factor .* df.baseline_utility
end


"Ensure the data/, tables/ and figs/ output directories exist under `rootdir`."
function ensure_dirs(rootdir)
    for dir in ("data", "tables", "figs")
        mkpath(joinpath(rootdir, dir))
    end
end


"Warm up Julia so experiment timings aren't polluted by first-call compilation."
function warmup()
    println("\nWARMING UP THE ENGINE")
    pop = generate_instance(10, 0:0.1:1, 1:10)
    greedy(pop; T=2, G=5)
    approximate(pop; T=2, G=5, K=15)
    println("READY TO RUMBLE\n")
end


## GLOBAL PARAMETERS

const UTIL_UPPER_BOUND = 50
const PILOT_BUDGETS = collect(2:4:30)      # {2, 6, ..., 30} for the pilot experiments
const SYNTHETIC_BUDGETS = collect(2:2:12)  # {2, 4, ..., 12} for the synthetic experiments


## EXPERIMENTS

"Write an experiment's results CSV and LaTeX summary table under `rootdir`."
function save_results(df, num, rootdir)
    CSV.write(joinpath(rootdir, "data", "exp$(num)-data.csv"), df)
    table = create_summary(df)
    open(joinpath(rootdir, "tables", "exp$(num)-summary.tex"), "w") do io
        show(io, "text/latex", table)
    end
end


"""
Build `reps` synthetic populations of size `n` by fitting a Normal distribution
to the pilot study's baseline utilities and drawing health probabilities from
`Uniform(0.5, 1)`.
"""
function synthetic_populations(rootdir, n, reps)
    _, baseline_utilities = extract_population(joinpath(rootdir, "pilotdata.csv"), UTIL_UPPER_BOUND)
    utility_distribution = fit(Normal, baseline_utilities)
    health_probs = Uniform(0.5, 1)
    return [generate_instance(n, health_probs, utility_distribution) for _ in 1:reps]
end

# Standard approx-vs-greedy algorithm pair, parameterised by the approximation
# parameter K.
approx_vs_greedy(K) = [
    (name=:approx, fn=approximate, args=Dict(:K => K, :verbose => false)),
    (name=:greedy, fn=greedy, args=Dict()),
]


"""
Run an approx-vs-greedy experiment on the pilot population (experiments 1–2).
Single population, so timings are kept clean by running single-threaded.
"""
function pilot_experiment(num; rootdir, G, budgets, K, seed, desc)
    println("\nSTARTING EXPERIMENT $(num): $(desc)")
    Random.seed!(seed)  # for reproducibility
    ensure_dirs(rootdir)
    trial_population, _ = extract_population(joinpath(rootdir, "pilotdata.csv"), UTIL_UPPER_BOUND)
    algs = approx_vs_greedy(K)
    df = run_experiments(algs, [trial_population], budgets, [G])
    add_comparisons!(df, algs)
    save_results(df, num, rootdir)
    return df
end


"""
Run an approx-vs-greedy experiment on `reps` synthetic populations of size `n`
(experiments 3–4). Run single-threaded: these use greedy's MOSEK conic
subroutine, and MOSEK's Linux aarch64 build aborts when called from a Julia
worker thread (it only runs safely on the main thread). Welfare/ratio plots are
written.
"""
function synthetic_experiment(num; rootdir, n, G, K, reps, seed, desc)
    println("\nSTARTING EXPERIMENT $(num): $(desc)")
    Random.seed!(seed)  # for reproducibility
    ensure_dirs(rootdir)
    populations = synthetic_populations(rootdir, n, reps)
    algs = approx_vs_greedy(K)
    df = run_experiments(algs, populations, SYNTHETIC_BUDGETS, [G])
    add_comparisons!(df, algs)
    save_results(df, num, rootdir)
    plot_welfares(df, joinpath(rootdir, "figs", "exp$(num)-welfares.pdf"))
    plot_ratios(df, joinpath(rootdir, "figs", "exp$(num)-ratios.pdf"))
    return df
end


# Experiment 1: approx vs greedy for pilot data and G=5.
experiment1(; rootdir=".") = pilot_experiment(1; rootdir, G=5, budgets=PILOT_BUDGETS, K=25, seed=1002,
    desc="pilot data and G=5")

# Experiment 2: approx vs greedy for pilot data and G=10.
experiment2(; rootdir=".") = pilot_experiment(2; rootdir, G=10, budgets=PILOT_BUDGETS, K=25, seed=1003,
    desc="pilot data and G=10")

# Experiment 3: approx vs greedy with synthetic populations, n=150, G=5.
experiment3(; rootdir=".") = synthetic_experiment(3; rootdir, n=150, G=5, K=20, reps=20, seed=2001,
    desc="synthetic populations, n=150, G=5")

# Experiment 4: approx vs greedy with synthetic populations, n=150, G=10.
experiment4(; rootdir=".") = synthetic_experiment(4; rootdir, n=150, G=10, K=20, reps=20, seed=2002,
    desc="synthetic populations, n=150, G=10")


"""
Experiment 5: non-overlapping (disjoint) vs 2-overlapping testing with n=10, G=10.
"""
function experiment5(; rootdir=".")
    println("\nSTARTING EXPERIMENT 5: disjoint vs 2-overlap")
    Random.seed!(3001)  # for reproducibility
    ensure_dirs(rootdir)
    n = 10
    reps = 20
    small_populations = [generate_instance(n, 0:0.1:1, 1:3) for _ in 1:reps]
    small_budgets = [2, 3, 4, 5]
    algs = [
        (name=:disjoint, fn=exact, args=Dict(:k => 1)),
        (name=:two_overlap, fn=exact, args=Dict(:k => 2)),
    ]
    G = n
    exp5 = run_experiments(algs, small_populations, small_budgets, [G]; multithread=true)
    # Reverse the alg order so diff/ratio are two_overlap vs disjoint (2-overlap ≥ disjoint, so diff ≥ 0).
    add_comparisons!(exp5, [algs[2], algs[1]])
    CSV.write(joinpath(rootdir, "data", "exp5-data.csv"), exp5)

    # Welfares (violin)
    @df exp5 violin(
        :budget,
        :disjoint_welfare,
        side=:left,
        linewidth=0,
        label="Disjoint",
        legend=false,
        xticks=small_budgets,
        xlabel="Test budget",
        ylabel="Welfare",
        size=(400,300))
    @df exp5 dotplot!(:budget, :disjoint_welfare, side=:left, marker=(:black,stroke(0)), label="")
    @df exp5 violin!(:budget, :two_overlap_welfare, side=:right, linewidth=0, label="2-Overlap")
    @df exp5 dotplot!(:budget, :two_overlap_welfare, side=:right, marker=(:black,stroke(0)), label="")
    Plots.pdf(joinpath(rootdir, "figs", "exp5-welfares.pdf"))

    # Ratios (violin)
    @df exp5 violin(
        :budget,
        :ratio,
        linewidth=0,
        label="",
        legend=false,
        xlabel="Test budget",
        ylabel="Welfare ratio",
        xticks=small_budgets,
        size=(400,300))
    @df exp5 dotplot!(:budget, :ratio, marker=(:black,stroke(0)), label="")
    Plots.pdf(joinpath(rootdir, "figs", "exp5-ratios.pdf"))

    # Disjoint vs 2-overlap scatter (alternative to the violin plot)
    plot_disjoint_vs_overlap(exp5, joinpath(rootdir, "figs", "exp5-scatter.pdf"))

    # Summary table
    output = combine(
        groupby(exp5, :budget),
        :disjoint_welfare => mean,
        :disjoint_time => roundmean => "disjoint time",
        :two_overlap_welfare => mean,
        :two_overlap_time => roundmean => "2-overlap time",
    )
    open(joinpath(rootdir, "tables", "exp5-summary.tex"), "w") do io; show(io, "text/latex", output); end
    return exp5
end


## COMMAND LINE INTERFACE

const EXPERIMENTS = Dict(
    "1" => experiment1,
    "2" => experiment2,
    "3" => experiment3,
    "4" => experiment4,
    "5" => experiment5,
)

# Sorted experiment numbers for help text and default "run all".
const EXPERIMENT_NUMBERS = sort(collect(keys(EXPERIMENTS)))

function parse_commandline(args)
    s = ArgParseSettings(description="Run pooled testing experiments")
    @add_arg_table! s begin
        "--experiments"
            help = "Comma-separated experiments to run, e.g. \"1,3,5\". Available: $(join(EXPERIMENT_NUMBERS, ", ")). Default: all."
            arg_type = String
            default = join(EXPERIMENT_NUMBERS, ",")
        "--rootdir"
            help = "Root directory for input data and output (default: current directory)."
            arg_type = String
            default = "."
        "--force"
            help = "Re-run experiments even if their output CSV already exists (default: skip completed experiments so interrupted runs can resume)."
            action = :store_true
    end
    return parse_args(args, s)
end

"Path to an experiment's output CSV; its existence marks the experiment complete."
output_csv(rootdir, key) = joinpath(rootdir, "data", "exp$(key)-data.csv")

function (@main)(args)
    parsed = parse_commandline(args)
    rootdir = parsed["rootdir"]
    force = parsed["force"]
    requested = split(parsed["experiments"], [',', ' ']; keepempty=false)

    warmup()

    for exp in requested
        key = strip(exp)
        if haskey(EXPERIMENTS, key)
            # Skip experiments that already produced output, so an interrupted
            # run (e.g. a reclaimed spot instance) resumes where it left off.
            if !force && isfile(output_csv(rootdir, key))
                @info "Skipping experiment $(key): output already exists at $(output_csv(rootdir, key)). Use --force to re-run."
                continue
            end
            try
                EXPERIMENTS[key](rootdir=rootdir)
                @info "Completed experiment $(key)"
            catch e
                @error "Failed to run experiment $(key)" exception=(e, catch_backtrace())
            end
        else
            @warn "Unknown experiment: $(key). Available: $(join(EXPERIMENT_NUMBERS, ", "))"
        end
    end

    @info "All requested experiments completed!"
end
