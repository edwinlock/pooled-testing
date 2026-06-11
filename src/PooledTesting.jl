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
using ProgressMeter

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

# Per-thread Gurobi environments (not safe to share across concurrent solves).
const GRB_ENVS = Vector{Union{Nothing, Gurobi.Env}}(nothing, Threads.nthreads())
function grb_env()
    tid = Threads.threadid()
    GRB_ENVS[tid] === nothing && (GRB_ENVS[tid] = Gurobi.Env())
    return GRB_ENVS[tid]
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

end # module
