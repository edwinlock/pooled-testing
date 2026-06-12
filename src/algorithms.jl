# Test-allocation algorithms. Each returns `(welfare, pools, guarantee, guarantee_post)`:
# `guarantee` is the a-priori additive approximation bound (0 for the exact methods)
# and `guarantee_post` the certified post-hoc bound `dual bound - exact welfare`
# read off the solved MILP (0 where it doesn't apply). `guarantee` is a
# reproducible design parameter; `guarantee_post` is per-instance and usually
# tighter, since the solver typically beats its gap tolerance.

"""
    greedy(population; T, G=5, verbose=false)

Greedily build a non-overlapping allocation of up to `T` tests of size ≤ `G`,
adding the (conic-exact) optimal single test at each step and removing the people
it covers. Returns `(welfare, pools, 0.0, 0.0)`.
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
    return welfare(pools, population), pools, 0.0, 0.0
end

"""
    exact(population; k=1, T, G, verbose=false)

Compute an optimal allocation of `T` tests of size ≤ `G` with each individual in
at most `k` tests (`k=1` is non-overlapping). Uses MOSEK for a single conic test,
else a Gurobi MILP. Returns `(welfare, pools, 0.0, 0.0)`.
"""
function exact(population::Population; k=1, T, G, verbose=false)
    pop = copy(population)
    remove_zeros!(pop)
    isempty(pop) && return 0.0, [], 0.0, 0.0   # (welfare, pools, guarantee, post) contract
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
    return welfare(pools, pop), pools, 0.0, 0.0
end

"Optimal single test via the conic (MOSEK) formulation."
conic(pop::Population; G=5, verbose=false) = exact(pop; T=1, G=G, verbose=verbose)

"""
    approximate(population; T, G=5, K=15, accuracy=nothing, start=nothing, verbose=false)

Near-optimal non-overlapping allocation via the clustered MILP, whose exponential
constraints are approximated by a `K`-segment piecewise-linear function. Returns
`(welfare, pools, guarantee, guarantee_post)`, where `guarantee` is the additive approximation
bound `T·ε` and `guarantee_post` the certified post-hoc bound `dual bound −
welfare` (per-instance, usually tighter — the solver typically beats its gap
tolerance). The returned solution is always the MILP's.

`start` is an optional candidate allocation (pools of person ids, e.g. greedy's
output). It is augmented to a complete MIP warm start (benchmarked ≈1.8x
faster) and, with `accuracy`, its welfare scales the error budget — any
feasible allocation is a valid lower bound on the optimum.

If `accuracy` is given (a relative error budget, e.g. `0.002` for 0.2%), `K` is
ignored: the budget is split between the PWL segment count and Gurobi's
absolute MIP gap via `accuracy_params`, and the returned `guarantee` includes
the gap term, so it remains the total additive bound. Requires `start`.
"""
function approximate(population::Population; T, G=5, K=15, accuracy=nothing, start=nothing, verbose=false)
    accuracy === nothing || start !== nothing ||
        throw(ArgumentError("`accuracy` needs a candidate solution to scale the error budget: pass `start` (e.g. greedy pools)."))
    pop = copy(population)
    remove_zeros!(pop)
    isempty(pop) && return 0.0, [], 0.0, 0.0
    clusters = pop2clusters(pop)
    q, u, n, keylist = cluster2vec(clusters)
    T = min(T, sum(n))                    # can't have more tests than people
    gapabs = 0.0
    if accuracy !== nothing
        K, gapabs = accuracy_params(q, u; T=T, G=G, target=accuracy*welfare(start, pop))
    end
    m, x, guarantee = approx_disjoint_model(q, u, n; T=T, G=G, K=K, verbose=verbose)
    start === nothing || full_clustered_start_vals(m, start, pop, keylist, q, u; G=G, K=K)
    if accuracy !== nothing
        set_optimizer_attribute(m, "MIPGapAbs", gapabs)
        set_optimizer_attribute(m, "MIPGap", 0.0)  # only the absolute criterion binds
        guarantee += gapabs
    end
    optimize!(m)
    _, p = retrieve(m, x, T, length(n))
    final_pools = uncluster(p, clusters, keylist)
    final_welfare = welfare(final_pools, pop)
    # Post-hoc certificate: the dual bound caps the model optimum, which caps the
    # true optimum (the PWL upper-approximates exp), and final_welfare is exact.
    guarantee_post = max(objective_bound(m) - final_welfare, 0.0)
    return final_welfare, final_pools, guarantee, guarantee_post
end
