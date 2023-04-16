using Combinatorics, ProgressBars, StaticArrays

include("iterator.jl")

"""
Compute all ways to form B disjoint pools from the set of pools given,
and return them as index lists.
"""
function disjoint_allocations(pools, B)
    B == 1 && return eachindex(pools)
    allocations = Vector{Int}[]
    # Iterate over all possible selections of B indices of pools
    for indices in combinations(eachindex(pools), B)
        selected_pools = pools[indices]  # recover pools associated with indices
        # Check selected pools are disjoint
        pairwise_disjoint = all(isdisjoint(pool, other) for (pool, other) in combinations(selected_pools, 2))
        pairwise_disjoint && push!(allocations, indices)
    end
    return allocations
end


"""
Given q and u, compute the utilities for all pools in `pools' and store them
in `pool_utilities`.
"""
@inline function compute_pool_utilities!(pool_utilities, pools, q, u)
    for (i, pool) in enumerate(pools)
        @inbounds pool_utilities[i] = prod(q[j] for j in pool) * sum(u[j] for j in pool)
    end
end


@inline function maximise_welfare(q, u, allocations, pool_utilities)
    max_welfare, welfare = -1.0, 0.0
    for a in allocations
        welfare = 0.0  # reset for each allocation
        for pool in a
            @inbounds welfare += pool_utilities[pool]
        end
        (welfare > max_welfare) && (max_welfare = welfare)
    end
    return max_welfare
end


function search(n::Int, probabilities, utilities)
    pools = collect(powerset(1:n, 1, n)) :: Vector{Vector{Int}}
    allocations1 = disjoint_allocations(pools, 1)
    allocations2 = disjoint_allocations(pools, 2)
    allocations3 = disjoint_allocations(pools, 3)
    # Preallocate vectors for computations in loops
    pool_utilities = similar(pools, Float64)  # create an array of Float64s with same length as pools
    max_diff, max_q, max_u = -1, fill(1.0,n), fill(1,n)
    prob_choices = collect(LexiIter(probabilities, n))
    @Threads.threads for q ∈ ProgressBar(prob_choices)
        for u ∈ Base.product(ntuple(_ -> utilities, n)...)
            # Compute utility achievable by each pool
            compute_pool_utilities!(pool_utilities, pools, q, u)
            # Compute maximum achievable welfare for budgets B ∈ {1,2,3}
            w1 = maximise_welfare(q, u, allocations1, pool_utilities)
            w2 = maximise_welfare(q, u, allocations2, pool_utilities)
            w3 = maximise_welfare(q, u, allocations3, pool_utilities)
            # Compute marginal return difference, and keep track of solution if it's the current maximum
            diff = w1 - 2*w2 + w3
            if diff > max_diff + 1e-10
                max_diff = diff
                max_q .= q
                max_u .= u
            end
        end
    end
end

probabilities = 0.1:0.1:1.0
utilities = 1:10
n = 6
@time search(n, probabilities, utilities)
