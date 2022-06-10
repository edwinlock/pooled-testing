using CSV, DataFrames, Statistics, Dates, StatsPlots
using ProgressBars

include("optimisation.jl")

"""
Run experiments with specified algorithms, populations, budgets and poolsizes.
"""
function run_experiments(algs, populations, budgets, poolsizes)
    df = DataFrame(:budget=>Int[], :population=>Int[], :poolsize => Int[])  # for storing results
    for alg in algs  # add columns for all algorithms
        df[!, "$(alg.name)_welfare"] = Float64[]
        df[!, "$(alg.name)_error"] = Float64[]
        df[!, "$(alg.name)_time"] = Millisecond[]
    end
    for (i, pop) in ProgressBar(enumerate(populations))
        for T in budgets
            # println("Testing budget: $(T).")
            for G in poolsizes
                result = Dict{Symbol, Any}(:budget => T, :population => i, :poolsize => G)
                for (name, fn, args) in algs
                    start = Dates.now()
                    w, pools, error = fn(pop; T=T, G=G, args...)
                    time = Dates.now() - start
                    result[Symbol("$(name)_welfare")] = w
                    result[Symbol("$(name)_error")] = error
                    result[Symbol("$(name)_time")] = time
                end
                push!(df, result)
            end
        end
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
roundmean(x) = sum(x) ÷ length(x)
function create_table(df; n=nothing, G=nothing)
    output = combine(
        groupby(df, :budget),
        :approx_welfare => mean,
        :approx_error => mean => "Add. error",
        :approx_time => roundmean => "approx time",
        :greedy_welfare => mean,
        :greedy_time => roundmean => "greedy time",
    )
    CSV.write("tables/summary-$(Dates.now())-$(n)-$(G).csv", output)
    return output
end

"""Plot welfares of approx and greedy algorithm from df."""
function plot_welfares(df; n=nothing, G=nothing, ylims=nothing)
    # Plot welfares of approx and greedy
    @df df violin(
        :budget,
        :approx_welfare,
        side=:left,
        linewidth=0,
        label="",
        legend=false,
        xticks=2:2:12,
        # title="Approx vs Greedy with n=$(n) and G=$(G).",
        xlabel="Test budget", ylabel="Welfare",
        ylims=ylims,
        size=(400,300))
    @df df dotplot!(:budget, :approx_welfare, side=:left, marker=(:black,stroke(0)), label="")
    @df df violin!(:budget, :greedy_welfare, side=:right, linewidth=0, label="Greedy")
    @df df dotplot!(:budget, :greedy_welfare, side=:right, marker=(:black,stroke(0)), label="")
    Plots.pdf("figs/welfare-$(Dates.now())-$(n)-$(G).pdf")
end


"""Plot welfare ratios of approx and greedy algorithm from df."""
function plot_ratios(df; n=nothing, G=nothing)
    @df df violin(
        :budget,
        :ratio,
        linewidth=0,
        label="",
        # label="Ratio",
        legend=:right,
        # title="Approx vs Greedy (ratio) with n=$(n) and G=$(G).",
        xlabel="Test budget",
        ylabel="Welfare ratio",
        xticks=2:2:12,
        ylims=ylims,
        size=(400,300)
        )
    @df df dotplot!(:budget, :ratio, marker=(:black,stroke(0)), label="")
    Plots.pdf("figs/ratios-$(Dates.now())-$(n)-$(G).pdf")
end


## BEGIN EXPERIMENTS
# Instructions: Comment out any experiments you don't want to run!

# Global constants for all experiments
health_probs = 0:0.1:1
utils = 1:10
budgets = [2, 4, 6, 8, 10, 12]
reps = 20
large_populations = [generate_instance(250, health_probs, utils) for _ in 1:reps]

# Dummy warmup for Julia to pre-compile all functions and get accurate running times for experiments
println("WARMING UP THE ENGINE")
pop = generate_instance(10, 0:0.1:1, 1:10)
greedy(pop; T=2, G=5)
approximate(pop; T=2, G=5, K=15)
println("READY TO RUMBLE")

# Experiment 1: Greedy vs Non-overlapping with n=250, G=5
println("STARTING EXPERIMENT 1")
n = 250
G = 5
algs = [(name=:approx, fn=approximate, args=Dict(:K => 18)), (name=:greedy, fn=greedy, args=Dict())]
exp1 = run_experiments(algs, large_populations, budgets, [G])
add_comparisons!(exp1, algs)
CSV.write("data/experiment1-$(Dates.now()).csv", exp1)
plot_welfares(exp1, n=n, G=G)
plot_ratios(exp1, n=n, G=G)
exp1_table = create_table(exp1, n=n, G=G)

# Experiment 2: Greedy vs Non-overlapping with n=250, G=10
println("STARTING EXPERIMENT 2")
n = 250
G = 10
algs = [(name=:approx, fn=approximate, args=Dict(:K => 21)), (name=:greedy, fn=greedy, args=Dict())]
exp2 = run_experiments(algs, large_populations, budgets, [G])
add_comparisons!(exp2, algs)
CSV.write("data/experiment2-$(Dates.now()).csv", exp2)
plot_welfares(exp2, n=n, G=G)
plot_ratios(exp2, n=n, G=G)
exp2_table = create_table(exp2, n=n, G=G)

# Experiment 3: Greedy vs Non-overlapping with n=250, G=250
println("STARTING EXPERIMENT 3")
n = 250
G = 250
algs = [(name=:approx, fn=approximate, args=Dict(:K => 65)), (name=:greedy, fn=greedy, args=Dict())]
exp3 = run_experiments(algs, large_populations, budgets, [G])
add_comparisons!(exp3, algs)
CSV.write("data/experiment3-$(Dates.now()).csv", exp3)
plot_welfares(exp3, n=n, G=G)
plot_ratios(exp3, n=n, G=G)
exp3_table = create_table(exp3, n=n, G=G)

# # Experiment 4: Non-overlapping vs 2-Overlapping with n=10, G=3
println("STARTING EXPERIMENT 4")
n = 10
reps = 20
small_populations = [generate_instance(n, 0:0.1:1, 1:3) for _ in 1:reps]
budgets = [2,3,4]
G = n
algs = [(name=:disjoint, fn=exact, args=Dict(:k => 1)), (name=:two_overlap, fn=exact, args=Dict(:k => 2))]
exp4 = run_experiments(algs, small_populations, budgets, [G])
add_comparisons!(exp4, [algs[2], algs[1]])

# Plot welfares for Experiment 4
@df exp4 violin(
    :budget,
    :disjoint_welfare,
    side=:left,
    linewidth=0,
    label="Disjoint",
    legend=false,
    xticks=budgets,
    xlabel="Test budget",
    ylabel="Welfare",
    # ylims=ylims,
    size=(400,300))
@df exp4 dotplot!(:budget, :disjoint_welfare, side=:left, marker=(:black,stroke(0)), label="")
@df exp4 violin!(:budget, :two_overlap_welfare, side=:right, linewidth=0, label="2-Overlap")
@df exp4 dotplot!(:budget, :two_overlap_welfare, side=:right, marker=(:black,stroke(0)), label="")
Plots.pdf("figs/welfare-small-$(Dates.now())-$(n)-$(G).pdf")

# Plot ratios for Experiment 4
@df exp4 violin(
    :budget,
    :ratio,
    linewidth=0,
    label="",
    legend=false,
    xlabel="Test budget",
    ylabel="Welfare ratio",
    xticks=budgets,
    # ylims=ylims,
    size=(400,300))
@df exp4 dotplot!(:budget, :ratio, marker=(:black,stroke(0)), label="")
Plots.pdf("figs/ratios-small-$(Dates.now())-$(n)-$(G).pdf")

# Create summary table for Experiment 4
output = combine(
    groupby(exp4, :budget),
    :disjoint_welfare => mean,
    :disjoint_time => roundmean => "disjoint time",
    :two_overlap_welfare => mean,
    :two_overlap_time => roundmean => "2-overlap time",
)
CSV.write("tables/summary-small-$(Dates.now())-$(n)-$(G).csv", output)

# End Experiment 4

## END EXPERIMENTS
