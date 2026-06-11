# Test-allocation algorithms. Each returns `(welfare, pools, error)`, where
# `error` is the MILP's additive approximation guarantee (0 for the exact methods).

"""
    greedy(population; T, G=5, verbose=false)

Greedily build a non-overlapping allocation of up to `T` tests of size ≤ `G`,
adding the (conic-exact) optimal single test at each step and removing the people
it covers. Returns `(welfare, pools, 0.0)`.
"""
function greedy(population::Population; T, G=5, verbose=false)
    pop = copy(population)  # defensive copy
    remove_zeros!(pop)
    welfares, pools = [], []
    for _ in 1:T
        isempty(pop) && break
        w, p = conic(pop; G=G, verbose=verbose)
        pool = p[1]                      # the single (=only) pool
        push!(welfares, w); push!(pools, pool)
        filter!(x -> !in(x.first, pool), pop)  # remove tested individuals
    end
    return welfare(pools, population), pools, 0.0
end

"""
    exact(population; k=1, T, G, verbose=false)

Compute an optimal allocation of `T` tests of size ≤ `G` with each individual in
at most `k` tests (`k=1` is non-overlapping). Uses MOSEK for a single conic test,
else a Gurobi MILP. Returns `(welfare, pools, 0.0)`.
"""
function exact(population::Population; k=1, T, G, verbose=false)
    pop = copy(population)
    remove_zeros!(pop)
    isempty(pop) && return 0.0, [], 0.0   # (welfare, pools, error) contract
    q, u, keylist = pop2vec(pop)
    m, x = if T == 1 && k == 1
        single_mosek(q, u; G=G, verbose=verbose)
    elseif k == 1
        milp_disjoint_model(q, u; T=T, G=G, verbose=verbose)
    else
        milp_overlap_model(q, u; k=k, T=T, G=G, verbose=verbose)
    end
    optimize!(m)
    _, p = retrieve(m, x, T, length(q))
    pools = [[keylist[i] for i in pool] for pool in p]
    return welfare(pools, pop), pools, 0.0
end

"Optimal single test via the conic (MOSEK) formulation."
conic(pop::Population; G=5, verbose=false) = exact(pop; T=1, G=G, verbose=verbose)

"""
    approximate(population; T, G=5, K=15, verbose=false)

Near-optimal non-overlapping allocation via the clustered MILP, whose exponential
constraints are approximated by a `K`-segment piecewise-linear function. Returns
`(welfare, pools, guarantee)`, where `guarantee` is the additive approximation
bound `T·ε`.
"""
function approximate(population::Population; T, G=5, K=15, verbose=false)
    pop = copy(population)
    remove_zeros!(pop)
    isempty(pop) && return 0.0, [], 0.0
    clusters = pop2clusters(pop)
    q, u, n, keylist = cluster2vec(clusters)
    T = min(T, sum(n))                    # can't have more tests than people
    m, x, guarantee = approx_disjoint_model(q, u, n; T=T, G=G, K=K, verbose=verbose)
    optimize!(m)
    _, p = retrieve(m, x, T, length(n))
    final_pools = uncluster(p, clusters, keylist)
    return welfare(final_pools, pop), final_pools, guarantee
end
