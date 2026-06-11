# Analyse the pooled-testing solve store (data/solves.db): build summary tables
# and plots from whatever solves are present. Safe to run while `run.jl` is still
# writing — it reads the store, it does not solve anything.
#
#     julia --project=. analyse.jl                  # summaries + plots for all experiments
#     julia --project=. analyse.jl --experiments 3  # a subset
#
# Output: data/expN-data.csv, tables/expN-summary.tex, figs/expN-*.pdf

using CSV, DataFrames, Statistics, Dates, StatsPlots, PrettyTables, ArgParse, SQLite, DBInterface, SHA, Pkg

const GUROBI_MIPGAP = 1e-4  # only needed so datastore.jl's param_key is consistent
include("datastore.jl")


## RESHAPING THE STORE

"""
Long-format solves for one experiment, pivoted to one row per (population, budget)
with a `<alg>_welfare` / `<alg>_time` column per algorithm — the shape the
summaries and plots expect. Adds `diff`/`ratio` between the two algorithms (the
first listed is the numerator).
"""
function experiment_frame(rootdir, experiment; algorder=nothing)
    df = filter(:experiment => ==(experiment), load_solves(rootdir))
    isempty(df) && return df
    algs = algorder === nothing ? sort(unique(df.alg)) : algorder
    keys = [:pop_index, :budget, :poolsize]
    wide  = unstack(df, keys, :alg, :welfare;   renamecols=a->Symbol("$(a)_welfare"))
    times = unstack(df, keys, :alg, :time_ms;   renamecols=a->Symbol("$(a)_time"))
    guars = unstack(df, keys, :alg, :guarantee; renamecols=a->Symbol("$(a)_guarantee"))
    wide = innerjoin(wide, times, guars, on=keys)
    a, b = algs[1], algs[2]
    wide.diff  = wide[!, Symbol("$(a)_welfare")] .- wide[!, Symbol("$(b)_welfare")]
    wide.ratio = wide[!, Symbol("$(a)_welfare")] ./ wide[!, Symbol("$(b)_welfare")]
    return sort!(wide, [:budget, :pop_index])
end


## SUMMARY TABLES

"Mean over populations of welfare, guarantee and (rounded) time per budget."
function summary_table(wide, algs)
    cols = Any[]
    for a in algs
        push!(cols, Symbol("$(a)_welfare") => mean => "$(a) welfare")
        gcol = Symbol("$(a)_guarantee")
        gcol in propertynames(wide) && push!(cols, gcol => mean => "$(a) guarantee")
        push!(cols, Symbol("$(a)_time") => (t -> round(Int, mean(t))) => "$(a) time (ms)")
    end
    return combine(groupby(wide, :budget), cols...)
end

"Write a summary table to data/ (CSV) and tables/ (LaTeX), and print it."
function write_summary(rootdir, num, wide, algs)
    table = summary_table(wide, algs)
    CSV.write(joinpath(rootdir, "data", "exp$(num)-data.csv"), wide)
    open(joinpath(rootdir, "tables", "exp$(num)-summary.tex"), "w") do io
        show(io, "text/latex", table)
    end
    println("\nExperiment $(num):"); pretty_table(table)
    return table
end


## PLOTS

"Side-by-side welfare violins for two algorithms across budgets."
function plot_welfares(wide, a, b, filename; budgets=2:2:12)
    @df wide violin(:budget, cols(Symbol("$(a)_welfare")); side=:left, linewidth=0, label="",
        legend=false, xticks=budgets, xlabel="Test budget", ylabel="Welfare", size=(400,300))
    @df wide dotplot!(:budget, cols(Symbol("$(a)_welfare")); side=:left, marker=(:black,stroke(0)), label="")
    @df wide violin!(:budget, cols(Symbol("$(b)_welfare")); side=:right, linewidth=0, label=string(b))
    @df wide dotplot!(:budget, cols(Symbol("$(b)_welfare")); side=:right, marker=(:black,stroke(0)), label="")
    Plots.pdf(filename)
end

"Welfare-ratio violins across budgets."
function plot_ratios(wide, filename; budgets=2:2:12)
    @df wide violin(:budget, :ratio; linewidth=0, label="", legend=false,
        xlabel="Test budget", ylabel="Welfare ratio", xticks=budgets, size=(400,300))
    @df wide dotplot!(:budget, :ratio; marker=(:black,stroke(0)), label="")
    Plots.pdf(filename)
end

"Scatter of one algorithm's welfare against another, coloured/shaped by budget, with y=x."
function plot_scatter(wide, a, b, filename)
    xa, yb = Symbol("$(a)_welfare"), Symbol("$(b)_welfare")
    @df wide scatter(cols(xa), cols(yb); group=:budget, markershape=:auto, alpha=0.7,
        xlabel="$(a) welfare", ylabel="$(b) welfare", legend=:topleft, legendtitle="Budget", size=(400,300))
    lo, hi = extrema([wide[!, xa]; wide[!, yb]])
    plot!([lo, hi], [lo, hi]; line=:dash, color=:gray, label="y=x")
    Plots.pdf(filename)
end


## DRIVING THE ANALYSIS PER EXPERIMENT

# (algorithms, whether to make welfare/ratio plots, whether to make a scatter)
const ANALYSES = Dict(
    "1" => (algs=["approx","greedy"], plots=false, scatter=false, budgets=2:4:30),
    "2" => (algs=["approx","greedy"], plots=false, scatter=false, budgets=2:4:30),
    "3" => (algs=["approx","greedy"], plots=true,  scatter=false, budgets=2:2:12),
    "4" => (algs=["approx","greedy"], plots=true,  scatter=false, budgets=2:2:12),
    "5" => (algs=["two_overlap","disjoint"], plots=true, scatter=true, budgets=2:5),
)

function analyse(rootdir, key)
    spec = ANALYSES[key]
    wide = experiment_frame(rootdir, parse(Int, key); algorder=spec.algs)
    if isempty(wide)
        @warn "No solves stored for experiment $(key) yet."
        return
    end
    for d in ("data", "tables", "figs"); mkpath(joinpath(rootdir, d)); end
    write_summary(rootdir, key, wide, spec.algs)
    a, b = spec.algs
    if spec.plots
        plot_welfares(wide, a, b, joinpath(rootdir, "figs", "exp$(key)-welfares.pdf"); budgets=spec.budgets)
        plot_ratios(wide, joinpath(rootdir, "figs", "exp$(key)-ratios.pdf"); budgets=spec.budgets)
    end
    spec.scatter && plot_scatter(wide, a, b, joinpath(rootdir, "figs", "exp$(key)-scatter.pdf"))
end


## COMMAND LINE INTERFACE

function parse_commandline(args)
    s = ArgParseSettings(description="Analyse the pooled-testing solve store.")
    @add_arg_table! s begin
        "--experiments"
            help = "Comma-separated experiments to analyse. Default: all."
            arg_type = String
            default = join(sort(collect(keys(ANALYSES))), ",")
        "--rootdir"
            arg_type = String
            default = "."
    end
    return parse_args(args, s)
end

function main(args)
    parsed = parse_commandline(args)
    rootdir = parsed["rootdir"]
    for key in strip.(split(parsed["experiments"], [',', ' ']; keepempty=false))
        haskey(ANALYSES, key) ? analyse(rootdir, key) : @warn("Unknown experiment: $(key)")
    end
end

# Run only when executed as a script (not when included, e.g. by the tests).
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
