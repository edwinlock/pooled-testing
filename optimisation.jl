using JuMP, Gurobi
using DataStructures

# Define constants
const GRB_ENV = Gurobi.Env()
const Population{T} = Dict{T, Tuple{Float64, Int}} where T  # interpretation: participant ids are mapped to (q, u)

include("utils.jl")
include("approximation-models.jl")
include("linear-models.jl")
include("conic-models.jl")


function greedy(pop::Population; T, G=5, verbose=false)
    pop = scale_utilities(pop, 15)  # makes defensive copy
    remove_zeros!(pop)
    isempty(pop) && return 0, Dict()
    welfares, pools, times = [], [], []
    start = 
    for t in 1:T
        w, p = conic(pop, G=G, verbose=verbose)
        pool = p[1]  # get the first (=only) pool
        push!(welfares, w)  # record welfare
        push!(pools, pool)  # record pool
        # Remove people in pool from population
        filter!(x->!in(x.first,pool), pop)
    end
    return sum(welfares), pools, 0.
end


function exact(pop::Population; k=1, T, G, verbose=false)
    pop = scale_utilities(pop, 15)  # makes defensive copy
    remove_zeros!(pop)
    isempty(pop) && return 0, Dict()
    q, u, keylist = pop2vec(pop)  # Get input vectors for model
    if T==1 && k==1
        m, x = single_mosek(q, u; G=G, verbose=verbose)
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

function approximate(pop::Population; T, G=5, K=15, verbose=false)
    pop = scale_utilities(pop, 15)  # makes defensive copy
    remove_zeros!(pop)
    isempty(pop) && return 0, Dict()
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
