using JuMP, Gurobi, DataStructures

const Population{T} = Dict{T, Tuple{Float64, Int}} where T  # interpretation: participant ids are mapped to (q, u)

# Apply Gurobi solver settings to a model. GUROBI_THREADS and GUROBI_MIPGAP are
# normally defined at the top of run.jl (before this file is included);
# fall back to sensible defaults if optimisation.jl is used standalone.
if !isdefined(@__MODULE__, :GUROBI_THREADS)
    const GUROBI_THREADS = 8
end
if !isdefined(@__MODULE__, :GUROBI_MIPGAP)
    # Loose convenience default for standalone use; run.jl/analyse.jl set their own
    # (1e-4). Comparison-grade code should set MIPGap explicitly, since a 0.1%
    # gap is too loose to resolve the MILP-vs-greedy difference (see log.md).
    const GUROBI_MIPGAP = 1e-3   # 0.1%
end
set_gurobi_params!(m) = set_optimizer_attributes(m,
    "Threads" => GUROBI_THREADS, "MIPGap" => GUROBI_MIPGAP)

# Gurobi environments are not safe to share across concurrent optimizations, so
# we keep one environment per thread and create them lazily on first use. Model
# code should call `grb_env()` instead of referencing a shared global.
const GRB_ENVS = Vector{Union{Nothing, Gurobi.Env}}(nothing, Threads.nthreads())
function grb_env()
    tid = Threads.threadid()
    if GRB_ENVS[tid] === nothing
        GRB_ENVS[tid] = Gurobi.Env()
    end
    return GRB_ENVS[tid]
end

include("models/utils.jl")
include("models/approximation-models.jl")
include("models/linear-models.jl")
include("models/conic-models.jl")


function greedy(population::Population; T, G=5, verbose=false)
    pop = copy(population)  # makes defensive copy
    remove_zeros!(pop)
    welfares, pools = [], []
    for t in 1:T
        isempty(pop) && break
        w, p = conic(pop, G=G, verbose=verbose)
        pool = p[1]  # get the first (=only) pool
        push!(welfares, w)  # record welfare
        push!(pools, pool)  # record pool
        # Remove people in pool from population
        filter!(x->!in(x.first,pool), pop)
    end
    w = welfare(pools, population)
    return w, pools, 0.
end


function exact(population::Population; k=1, T, G, verbose=false)
    pop = copy(population)  # makes defensive copy
    remove_zeros!(pop)
    isempty(pop) && return 0., [], 0.  # (welfare, pools, error) contract
    q, u, keylist = pop2vec(pop)  # Get input vectors for model
    if T==1 && k==1
        m, x = single_mosek(q, u; G=G, verbose=verbose)
    elseif k==1
        m, x = milp_disjoint_model(q, u; T=T, G=G, verbose=verbose)
    else
        m, x = milp_overlap_model(q, u; k=k, T=T, G=G, verbose=verbose)
    end
    optimize!(m)
    _, p = retrieve(m, x, T, length(q))
    pools = [[keylist[i] for i in pool] for pool in p]
    w = welfare(pools, pop)
    return w, pools, 0.
end


conic(pop::Population; G=5, verbose=false) = exact(pop, T=1, G=G, verbose=verbose)


function approximate(population::Population; T, G=5, K=15, verbose=false)
    pop = copy(population)  # makes defensive copy
    remove_zeros!(pop)
    isempty(pop) && return 0., [], 0.  # (welfare, pools, error) contract
    clusters = pop2clusters(pop)
    q, u, n, keylist = cluster2vec(clusters)
    T = min(T,sum(n))  # can't have more tests than people
    m, x, error = approx_disjoint_model(q, u, n; T=T, G=G, K=K, verbose=verbose)
    optimize!(m)
    _, p = retrieve(m, x, T, length(n))
    final_pools = uncluster(p, clusters, keylist)
    w = welfare(final_pools, pop)
    return w, final_pools, error
end


# function overlap(pop::Population; k, T, G=5, K=15, verbose=false)
#     pop = scale_utilities(pop, 15)  # makes defensive copy
#     remove_zeros!(pop)
#     isempty(pop) && return 0, Dict()
#     q, u, keylist = pop2vec(pop)  # Get input vectors for model
#     n = length(q)
#     T = min(T,n)  # can't have more tests than people
#     m, x, error = approx_model(q, u; k=k, T=T, G=G, K=K, verbose=verbose)
#     optimize!(m)
#     _, p = retrieve(m, x, T, length(q))
#     pools = [[keylist[i] for i in pool] for pool in p]
#     w = welfare(pools, pop)
#     return w, pools, error
# end
