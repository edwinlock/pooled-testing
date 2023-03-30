# Install all packages that aren't yet installed
# using Pkg
# Pkg.add([
#     "CSV", "DataFrames", "Statistics", "Dates", "StatsPlots",
#     "Distributions", "ProgressBars", "DataStructures", "Gurobi"])

# Load packages
using CSV, DataFrames, Statistics, Dates, StatsPlots, Distributions, ProgressBars, DataStructures, Gurobi

# Include optimisation code
include("optimisation.jl")

"""
Convenience function to run experiments with specified algorithms, populations, budgets and poolsizes.
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
function create_summary(df)
    output = combine(
        groupby(df, :budget),
        :approx_welfare => mean,
        :approx_error => mean => "Add. error",
        :approx_time => roundmean => "approx time",
        :greedy_welfare => mean,
        :greedy_time => roundmean => "greedy time",
    )
    # CSV.write("tables/summary-$(Dates.now())-$(n)-$(G).csv", output)
    return output
end


"""Plot welfares of approx and greedy algorithm from df."""
function plot_welfares(df, filename)
    # Plot welfares of approx and greedy
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


"Extract population data from CSV file"
function extract_population(filename, utility_upper_bound)
    df = DataFrame(CSV.File(filename))
    scaling_factor = utility_upper_bound / maximum(df.baseline_utility)
    trial_population = Population{Int}()  # with integral utilities 
    for row in eachrow(df)
        q = row.health_probability
        u = Int(round(row.baseline_utility*scaling_factor))
        trial_population[row.id] = (q, u)
    end
    return trial_population, scaling_factor .* df.baseline_utility
end


# Dummy warmup for Julia to pre-compile all functions and get accurate running times for experiments
println("\nWARMING UP THE ENGINE")
pop = generate_instance(10, 0:0.1:1, 1:10)
greedy(pop; T=2, G=5)
approximate(pop; T=2, G=5, K=15)
println("READY TO RUMBLE\n")


## BEGIN EXPERIMENTS
# Instructions: Comment out any experiments you don't want to run!

# Set global parameters
util_upper_bound = 50
budgets = [2, 4, 6, 8, 10, 12]

# Extract pilot data from CSV file
trial_population, baseline_utilities = extract_population("pilotdata.csv", util_upper_bound)
n = length(trial_population)

# Set budget and algorithm parameters for experiments
algs = [
    (name=:approx, fn=approximate, args=Dict(:K => 25, :verbose => false)),
    (name=:greedy, fn=greedy, args=Dict())
]

# Experiment 0 for Evi: Greedy vs Non-overlapping for pilot data, G=5, and a testing budget of 30
# println("\nSTARTING EXPERIMENT 0 for Evi: pilot data, G=5, and a testing budget of 30.")
# G = 5
# exp0 = run_experiments(algs, [trial_population], [30], [G])  # note that I've set the budget to 30 here
# add_comparisons!(exp0, algs)
# CSV.write("data/exp0-data.csv", exp0)
# exp1_table = create_summary(exp0)
# open("tables/exp0-summary.tex", "w") do io; show(io, "text/latex", exp0_table); end

# Experiment 1: Greedy vs Non-overlapping for pilot data and G=5
println("\nSTARTING EXPERIMENT 1: pilot data and G=5")
G = 5
exp1 = run_experiments(algs, [trial_population], budgets, [G])
add_comparisons!(exp1, algs)
CSV.write("data/exp1-data.csv", exp1)
exp1_table = create_summary(exp1)
open("tables/exp1-summary.tex", "w") do io; show(io, "text/latex", exp1_table); end


# Experiment 2: Greedy vs Non-overlapping for pilot data and G=10
println("\nSTARTING EXPERIMENT 2: pilot data and G=10")
G = 10
exp2 = run_experiments(algs, [trial_population], budgets, [G])
add_comparisons!(exp2, algs)
CSV.write("data/exp2-data.csv", exp2)
exp2_table = create_summary(exp2)
open("tables/exp2-summary.tex", "w") do io; show(io, "text/latex", exp2_table); end

# # SYNTHETIC EXPERIMENTS
# # Fit distribution to utilities from pilot study
# utility_distribution = fit(Normal, baseline_utilities)
# # Verify that the distribution is a good fit:
# # histogram(baseline_utilities, normalize=true, bins=20, fill=:lightgray, fillalpha=0.3, legend=false)
# # density!(baseline_utilities, linecolor=:black)
# # utility_pdf(x) = Distributions.pdf(utility_distribution, x)
# # plot!(utility_pdf, linecolor=:blue, linestyle=:dash)
# n = 200  # population size 
# health_probs = Uniform(0.5,1)  # range of health probabilities to draw from
# reps = 20  # number of populations to generate
# populations = [generate_instance(n, health_probs, utility_distribution) for _ in 1:reps]
# budgets = [2, 4, 6, 8, 10, 12]

# # Experiment 3: Greedy vs Non-overlapping with population size n and G=5
# println("\nSTARTING EXPERIMENT 3")
# G = 5
# algs = [
#     (name=:approx, fn=approximate, args=Dict(:K => 20)),
#     (name=:greedy, fn=greedy, args=Dict())
# ]
# exp3 = run_experiments(algs, populations, budgets, [G])
# add_comparisons!(exp3, algs)
# CSV.write("data/exp3-data.csv", exp3)
# plot_welfares(exp3, "figs/exp3-welfares.pdf")
# plot_ratios(exp3, "figs/exp3-ratios.pdf")
# exp3_table = create_summary(exp3)
# open("tables/exp3-summary.tex", "w") do io; show(io, "text/latex", exp3_table); end

# # Experiment 4: Greedy vs Non-overlapping with population size n and G=10
# println("\nSTARTING EXPERIMENT 4")
# G = 10
# algs = [
#     (name=:approx, fn=approximate, args=Dict(:K =>20)),
#     (name=:greedy, fn=greedy, args=Dict())
# ]
# exp4 = run_experiments(algs, populations, budgets, [G])
# add_comparisons!(exp4, algs)
# CSV.write("data/exp4-data.csv", exp4)
# plot_welfares(exp4, "figs/exp4-welfares.pdf")
# plot_ratios(exp4, "figs/exp4-ratios.pdf")
# exp4_table = create_summary(exp4)
# open("tables/exp4-summary.tex", "w") do io; show(io, "text/latex", exp4_table); end


# # # Experiment 5: Greedy vs Non-overlapping with n=150, G=150
# # println("\nSTARTING EXPERIMENT 5")
# # G = 150
# # algs = [(name=:approx, fn=approximate, args=Dict(:K =>10)), (name=:greedy, fn=greedy, args=Dict())]
# # exp5 = run_experiments(algs, populations, budgets, [G])
# # add_comparisons!(exp5, algs)
# # CSV.write("data/experiment5-$(Dates.now()).csv", exp5)
# # plot_welfares(exp5, n=n, G=G)
# # plot_ratios(exp5, n=n, G=G)
# # exp5_table = create_summary(exp5, n=n, G=G)


# # ## SECTION: TOWARDS OVERLAPPING TESTING

# # # Experiment 6: Non-overlapping vs 2-Overlapping with n=10, G=3
# # println("\nSTARTING EXPERIMENT 6")
# # n = 10
# # reps = 20
# # small_populations = [generate_instance(n, 0:0.1:1, 1:3) for _ in 1:reps]
# # budgets = [2,3,4]
# # G = n
# # algs = [
# #     (name=:disjoint, fn=exact, args=Dict(:k => 1)),
# #     (name=:two_overlap, fn=exact, args=Dict(:k => 2))
# # ]
# # exp6 = run_experiments(algs, small_populations, budgets, [G])
# # add_comparisons!(exp6, [algs[2], algs[1]])
# # CSV.write("data/exp6-data.csv", exp6)

# # # Plot welfares for Experiment 6
# # @df exp6 violin(
# #     :budget,
# #     :disjoint_welfare,
# #     side=:left,
# #     linewidth=0,
# #     label="Disjoint",
# #     legend=false,
# #     xticks=budgets,
# #     xlabel="Test budget",
# #     ylabel="Welfare",
# #     # ylims=ylims,
# #     size=(400,300))
# # @df exp6 dotplot!(:budget, :disjoint_welfare, side=:left, marker=(:black,stroke(0)), label="")
# # @df exp6 violin!(:budget, :two_overlap_welfare, side=:right, linewidth=0, label="2-Overlap")
# # @df exp6 dotplot!(:budget, :two_overlap_welfare, side=:right, marker=(:black,stroke(0)), label="")
# # Plots.pdf("figs/exp6-welfares.pdf")

# # # Plot ratios for Experiment 6
# # @df exp6 violin(
# #     :budget,
# #     :ratio,
# #     linewidth=0,
# #     label="",
# #     legend=false,
# #     xlabel="Test budget",
# #     ylabel="Welfare ratio",
# #     xticks=budgets,
# #     # ylims=ylims,
# #     size=(400,300))
# # @df exp6 dotplot!(:budget, :ratio, marker=(:black,stroke(0)), label="")
# # Plots.pdf("figs/exp6-ratios.pdf")

# # # Create summary table for Experiment 6
# # output = combine(
# #     groupby(exp6, :budget),
# #     :disjoint_welfare => mean,
# #     :disjoint_time => roundmean => "disjoint time",
# #     :two_overlap_welfare => mean,
# #     :two_overlap_time => roundmean => "2-overlap time",
# # )
# # open("tables/exp6-summary.tex", "w") do io; show(io, "text/latex", output); end

# # # End Experiment 6

# ## END EXPERIMENTS
