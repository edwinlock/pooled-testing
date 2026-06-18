# Analyse the pooled-testing solve store (data/solves.db): build summary tables
# and plots from whatever solves are present. Safe to run while `run.jl` is still
# writing — it reads the store, it does not solve anything.
#
#     julia --project=. analyse.jl                  # summaries + plots for all experiments
#     julia --project=. analyse.jl --experiments 3  # a subset
#
# Output: data/expN-data.csv, tables/expN-summary.tex, figs/expN-*.pdf

using PooledTesting, CSV, DataFrames, Statistics, StatsPlots, PrettyTables, ArgParse, Printf

const ROOT = @__DIR__   # the data/ store and outputs live alongside this script

# Same experiment parameters as run.jl (honours POOLED_CONSTANTS).
include(abspath(get(ENV, "POOLED_CONSTANTS", joinpath(@__DIR__, "constants.jl"))))


## RESHAPING THE STORE

function select_params(df, algs)
    selected = Dict{String,String}()
    for alg in algs
        rows = filter(:alg => ==(alg), df)
        isempty(rows) && continue
        sort!(rows, :timestamp)
        selected[alg] = rows.params[end]
        if length(unique(rows.params)) > 1
            @info "Multiple parameter regimes for $(alg); using latest." params=selected[alg]
        end
    end
    return selected
end

function filter_params(df, selected)
    isempty(selected) && return df
    return filter(row -> haskey(selected, row.alg) && row.params == selected[row.alg], df)
end

"""
Long-format solves for one experiment, pivoted to one row per (population, budget)
with a `<alg>_welfare` / `<alg>_time` column per algorithm — the shape the
summaries and plots expect. Adds `diff`/`ratio` between the two algorithms (the
first listed is the numerator).
"""
function experiment_frame(rootdir, experiment; algorder=nothing, budgets=nothing)
    df = filter(:experiment => ==(experiment), load_solves(rootdir))
    isempty(df) && return df
    # Restrict to the budgets defined in constants.jl (via the experiment spec).
    # The store may hold extra budgets from earlier runs; analyse only the
    # configured ones, leaving the surplus rows in the store untouched.
    if budgets !== nothing
        df = filter(:budget => in(Set(budgets)), df)
        isempty(df) && return df
    end
    algs = algorder === nothing ? sort(unique(df.alg)) : algorder
    params = select_params(df, algs)
    df = filter_params(df, params)
    isempty(df) && return df
    keys = :n in propertynames(df) ? [:pop_index, :budget, :poolsize, :n] : [:pop_index, :budget, :poolsize]
    wide  = unstack(df, keys, :alg, :welfare;        renamecols=a->Symbol("$(a)_welfare"))
    times = unstack(df, keys, :alg, :time_ms;        renamecols=a->Symbol("$(a)_time"))
    guars = unstack(df, keys, :alg, :guarantee;      renamecols=a->Symbol("$(a)_guarantee"))
    posts = unstack(df, keys, :alg, :guarantee_post; renamecols=a->Symbol("$(a)_guarantee_post"))
    wide = innerjoin(wide, times, guars, posts, on=keys)
    for alg in algs
        haskey(params, alg) && (wide[!, Symbol("$(alg)_params")] = fill(params[alg], nrow(wide)))
    end
    a, b = algs[1], algs[2]
    wide.diff  = wide[!, Symbol("$(a)_welfare")] .- wide[!, Symbol("$(b)_welfare")]
    wide.ratio = wide[!, Symbol("$(a)_welfare")] ./ wide[!, Symbol("$(b)_welfare")]
    return sort!(wide, [:budget, :pop_index])
end


## SUMMARY TABLES

# Mean ignoring missing cells (a live run hasn't filled every solve yet);
# `missing` when a whole group is missing.
meanskip(c) = (v = collect(skipmissing(c)); isempty(v) ? missing : mean(v))
intskip(c) = (v = meanskip(c); v === missing ? missing : round(Int, v))

"Mean over populations of welfare, guarantees and (rounded) time per budget."
function summary_table(wide, algs)
    cols = Any[]
    for a in algs
        push!(cols, Symbol("$(a)_welfare") => meanskip => "$(a) welfare")
        for (suffix, label) in (("guarantee", "guarantee"), ("guarantee_post", "post-hoc"))
            gcol = Symbol("$(a)_$(suffix)")
            if gcol in propertynames(wide) && a != "greedy"
                push!(cols, gcol => meanskip => "$(a) $(label)")
            end
        end
        push!(cols, Symbol("$(a)_time") => intskip => "$(a) time (ms)")
    end
    :diff in propertynames(wide) && push!(cols, :diff => meanskip => "welfare diff")
    :ratio in propertynames(wide) && push!(cols, :ratio => meanskip => "welfare ratio")
    return combine(groupby(wide, :budget), cols...)
end

comma_int(x) = replace(string(round(Int, x)), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
fmt2(x) = x === missing ? "--" : comma_int(floor(Int, x)) * @sprintf("%.2f", x - floor(x))[2:end]
fmt6(x) = x === missing ? "--" : @sprintf("%.6f", x)

function trim_decimal(s)
    s = replace(s, r"0+$" => "")
    return replace(s, r"\.$" => "")
end

function human_time_ms(x)
    x === missing && return "--"
    x < 1_000 && return "$(comma_int(x)) ms"
    x < 60_000 && return "$(trim_decimal(@sprintf("%.2f", x / 1_000))) s"
    x < 3_600_000 && return "$(trim_decimal(@sprintf("%.2f", x / 60_000))) min"
    return "$(trim_decimal(@sprintf("%.2f", x / 3_600_000))) h"
end

function paper_exp1_table(wide)
    table = combine(groupby(wide, :budget),
        :approx_welfare => meanskip => :approx_welfare,
        :approx_guarantee => meanskip => :approx_guarantee,
        :approx_guarantee_post => meanskip => :approx_guarantee_post,
        :approx_time => meanskip => :approx_time,
        :greedy_welfare => meanskip => :greedy_welfare,
        :greedy_time => meanskip => :greedy_time)
    # Upper bound on optimal/greedy welfare, from the tighter (post-hoc) certificate.
    table.apx_to_optimal = (table.approx_welfare .+ table.approx_guarantee_post) ./ table.greedy_welfare
    return sort!(table, :budget)
end

function write_paper_exp1_tex(path, wide)
    table = paper_exp1_table(wide)
    G = first(skipmissing(wide.poolsize))
    n = :n in propertynames(wide) ? first(skipmissing(wide.n)) : 130
    budgets = join(table.budget, ", ")
    open(path, "w") do io
        println(io, raw"\begin{table}[tb!]")
        println(io, raw"    \centering")
        println(io, raw"    \begin{threeparttable}")
        println(io, "    \\caption{Performance of the MILP and \\greedy{} on pilot data (\$G=$(G)\$).}")
        println(io, raw"    \label{table:experiment1}")
        println(io, raw"    {\small")
        println(io, raw"    \renewcommand{\arraystretch}{1.15}")
        println(io, raw"    \begin{tabular}{@{} crccrrcr @{}}")
        println(io, raw"        \toprule")
        println(io, "        & \\multicolumn{4}{c}{{MILP}} & \\multicolumn{3}{c}{{Greedy}} \\\\")
        println(io, raw"        \cmidrule(lr){2-5}")
        println(io, raw"        \cmidrule(l){6-8}")
        println(io, "        {Budget} & Welfare & Guarantee & Post-hoc & Time & Welfare & Apx To Optimal & Time\\\\")
        println(io, raw"        \midrule")
        for row in eachrow(table)
            println(io, "        $(row.budget) & $(fmt2(row.approx_welfare)) & $(fmt2(row.approx_guarantee)) & $(fmt2(row.approx_guarantee_post)) & $(human_time_ms(row.approx_time)) & $(fmt2(row.greedy_welfare)) & $(fmt6(row.apx_to_optimal)) & $(human_time_ms(row.greedy_time)) \\\\")
        end
        println(io, raw"        \bottomrule")
        println(io, raw"    \end{tabular}")
        println(io, raw"    }")
        println(io, raw"    \begin{tablenotes}[flushleft]")
        println(io, raw"    \footnotesize")
        println(io, "    \\item[] \\emph{Note.} Summary showing welfare and computation time for the MILP and \\greedy{} on the pilot data with a population of \$n=$(n)\$ and pool size constraint \$G=$(G)\$, with testing budgets \$B \\in \\{$(budgets)\\}\$. Welfare figures are deterministic values computed on the pilot population. The column ``Guarantee'' reports the a-priori additive approximation guarantee of the MILP relative to optimal non-overlapping welfare; ``Post-hoc'' reports the certified per-instance bound computed after the solve as the solver's dual bound minus the exact MILP welfare, which is at most the guarantee. The column ``Apx To Optimal'' reports the upper bound on the ratio between optimal non-overlapping welfare and \\greedy{} welfare, computed as the MILP welfare plus the post-hoc bound, divided by \\greedy{} welfare.")
        println(io, raw"    \end{tablenotes}")
        println(io, raw"    \end{threeparttable}")
        println(io, raw"\end{table}")
    end
    return table
end

"Write a summary table to data/ (CSV) and tables/ (LaTeX), and print it."
function write_summary(rootdir, num, wide, algs)
    table = summary_table(wide, algs)
    CSV.write(joinpath(rootdir, "data", "exp$(num)-data.csv"), wide)
    tex_path = joinpath(rootdir, "tables", "exp$(num)-summary.tex")
    if string(num) == "1" && algs == ["approx", "greedy"]
        write_paper_exp1_tex(tex_path, wide)
    else
        open(tex_path, "w") do io
            show(io, "text/latex", table)
        end
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

"""
Per-budget histograms of the percentage welfare improvement of `a` over `b`
(`100*(ratio-1)`) across populations: one stacked panel per budget, sharing the
x-axis, with a dashed line at 0% (the two algorithms tie). Surfaces the full
spread the welfare violins hide — how many populations gain nothing vs a lot.
"""
function plot_gain_hist(wide, a, b, filename; budgets=nothing)
    bs = budgets === nothing ? sort(unique(wide.budget)) : sort(filter(in(Set(wide.budget)), budgets))
    gain = 100 .* (collect(skipmissing(wide.ratio)) .- 1)
    lo, hi = extrema(gain)                        # shared x-range across panels
    edges = range(lo, hi; length=21)
    panels = map(enumerate(bs)) do (i, bud)
        g = 100 .* (collect(skipmissing(filter(:budget => ==(bud), wide).ratio)) .- 1)
        p = histogram(g; bins=edges, xlims=(lo, hi), legend=false, color=:steelblue,
            linecolor=:white, ylabel="B=$(bud)", yticks=nothing,
            xlabel=(i == length(bs) ? "$(a) over $(b) welfare (%)" : ""))
        vline!(p, [0.0]; line=(:dash, :gray), label="")
        p
    end
    plot(panels...; layout=(length(panels), 1), size=(420, 130*length(panels)),
        plot_title="Per-population welfare gain", link=:x)
    Plots.pdf(filename)
end


## DRIVING THE ANALYSIS PER EXPERIMENT

# (algorithms, whether to make welfare/ratio plots, whether to make a scatter).
# Budgets (used for plot xticks) come from the shared experiment constants.
const ANALYSES = Dict(
    "1" => (algs=["approx","greedy"], plots=false, scatter=false, budgets=PILOT_BUDGETS),
    "2" => (algs=["approx","greedy"], plots=false, scatter=false, budgets=PILOT_BUDGETS),
    "3" => (algs=["approx","greedy"], plots=true,  scatter=true, budgets=SYNTHETIC_BUDGETS),
    "4" => (algs=["approx","greedy"], plots=true,  scatter=true, budgets=SYNTHETIC_BUDGETS),
    "5" => (algs=["two_overlap","disjoint"], plots=true, scatter=true, budgets=[2,3,4,5]),
)

function analyse(rootdir, key)
    spec = ANALYSES[key]
    wide = experiment_frame(rootdir, parse(Int, key); algorder=spec.algs, budgets=spec.budgets)
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
        plot_gain_hist(wide, a, b, joinpath(rootdir, "figs", "exp$(key)-gainhist.pdf"); budgets=spec.budgets)
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
            default = ROOT
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
