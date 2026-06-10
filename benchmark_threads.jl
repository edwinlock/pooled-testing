# Benchmark Gurobi solve time vs thread count on a single fixed MILP.
#
# Solves the pilot-data MILP at one budget (default B=18, G=5, K=25) for a range
# of Gurobi thread counts and prints the wall-clock for each. Use it to find a
# good thread setting for this hardware (MILP branch-and-bound parallelises
# poorly, so more threads can be slower).
#
#   julia --project=. benchmark_threads.jl
#   julia --project=. benchmark_threads.jl 22          # different budget
#   julia --project=. benchmark_threads.jl 18 1,8,32   # custom thread list
#
# Note: does NOT need -t auto; this benchmarks Gurobi's *internal* threads, not
# Julia threads.

using CSV, DataFrames, JuMP, Gurobi, DataStructures, Dates, PrettyTables

# Pull in the model code (Population type, approx_disjoint_model, helpers).
# GUROBI_THREADS/MIPGAP get default values from optimisation.jl's fallback; we
# override Threads per solve below, and keep MIPGap at the paper's default.
include("optimisation.jl")

const UTIL_UPPER_BOUND = 50
const MIPGAP = 1e-4   # paper default; set explicitly so the benchmark never
                      # depends on optimisation.jl's fallback value

"Load the pilot population with integral utilities (mirrors experiments.jl)."
function load_pilot(filename)
    df = DataFrame(CSV.File(filename))
    scaling_factor = UTIL_UPPER_BOUND / maximum(df.baseline_utility)
    pop = Population{Int}()
    for row in eachrow(df)
        u = clamp(round(Int, row.baseline_utility * scaling_factor), 1, UTIL_UPPER_BOUND)
        pop[row.id] = (row.health_probability, u)
    end
    return pop
end

"Build the approximate MILP for `pop` at budget T, set `threads`, solve, return seconds."
function timed_solve(pop, T, G, K, threads)
    p = copy(pop); remove_zeros!(p)
    clusters = pop2clusters(p)
    q, u, n, _ = cluster2vec(clusters)
    T = min(T, sum(n))
    m, _, _ = approx_disjoint_model(q, u, n; T=T, G=G, K=K, verbose=false)
    # approx_disjoint_model already called set_gurobi_params! (Threads=GUROBI_THREADS,
    # MIPGap=fallback). Override both here, then read them back to confirm what
    # Gurobi will actually use.
    set_optimizer_attribute(m, "Threads", threads)
    set_optimizer_attribute(m, "MIPGap", MIPGAP)
    actual_threads = get_optimizer_attribute(m, "Threads")
    actual_gap = get_optimizer_attribute(m, "MIPGap")  # confirm the override took
    start = now()
    optimize!(m)
    elapsed = now() - start
    return elapsed, objective_value(m), MOI.get(m, MOI.RelativeGap()), actual_threads, actual_gap
end

# --- Parse args ---
budget      = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 14
threadlist  = length(ARGS) >= 2 ? parse.(Int, split(ARGS[2], ",")) : [1, 2, 4, 8, 16, 32, 64, 96]
G, K = 5, 25

pop = load_pilot("pilotdata.csv")
println("Benchmarking pilot MILP at B=$(budget), G=$(G), K=$(K), MIPGap=$(MIPGAP)")
println("Population size: $(length(pop))  |  machine cores: $(Sys.CPU_THREADS)\n")

# Warm up (compile + Gurobi init) on a trivial solve so timings are clean.
print("Warming up... "); timed_solve(pop, 2, G, K, 1); println("done")

# Confirm the actual MIPGap Gurobi will use (the build log can show a stale value).
_, _, _, _, actual_gap = timed_solve(pop, 2, G, K, 1)
println("Confirmed MIPGap in use: $(actual_gap)\n")

results = DataFrame(Threads=Int[], Used=Int[], var"Time (s)"=Float64[],
                    Objective=Float64[], RelGap=Float64[])
for t in threadlist
    elapsed, obj, gap, used, _ = timed_solve(pop, budget, G, K, t)
    secs = Dates.value(elapsed) / 1000  # ms -> seconds
    push!(results, (t, round(Int, used), round(secs; digits=1),
                    round(obj; digits=2), round(gap; sigdigits=3)))
end
pretty_table(results)
