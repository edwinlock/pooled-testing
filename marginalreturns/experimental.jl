using Combinatorics, ProgressBars

include("iterator.jl")


function allocations(n::Int)
    pools = collect(powerset(1:n, 1, n)) :: Vector{Vector{Int}}
    singles :: Vector{Int} = eachindex(pools)
    pairs = Tuple{Int,Int}[(i,j) for (i,j) in combinations(singles, 2) if isdisjoint(pools[i], pools[j])]
    triples = Tuple{Int,Int}[]
    for indices in combinations(eachindex(pools), 3)
        i, j, k = indices
        selected_pools = pools[indices]  # recover pools associated with indices
        # Check selected pools are disjoint
        pairwise_disjoint = all(isdisjoint(pool, other) for (pool, other) in combinations(selected_pools, 2))
        if pairwise_disjoint
            pair_index = findfirst(x->x==(i,j), pairs)  # find index of pair (i,j) in `pairs`
            push!(triples, (pair_index, k))
        end
    end
    return pools, singles, pairs, triples
end


function compute_utilities!(
    single_utils::Vector{Float64},
    pair_utils::Vector{Float64},
    triple_utils::Vector{Float64},
    pools::Vector{Vector{Int}},
    pairs::Vector{Tuple{Int,Int}},
    triples::Vector{Tuple{Int,Int}},
    q,
    u,
    )
    for (i, pool) in enumerate(pools)
        @inbounds single_utils[i] = prod(q[j] for j in pool) * sum(u[j] for j in pool)
    end
    for (i, (j,k)) in enumerate(pairs)
        @inbounds pair_utils[i] = single_utils[j] + single_utils[k]
    end
    for (i, (j,k)) in enumerate(triples)
        @inbounds triple_utils[i] = pair_utils[j] + single_utils[k]
    end
    return nothing
end


function faster_search(n::Int, probabilities, utilities)
    pools, singles, pairs, triples = allocations(n)
    # Preallocate vectors for computations in loops
    single_utils = similar(singles, Float64)
    pair_utils = similar(pairs, Float64)
    triple_utils = similar(triples, Float64)
    max_diff, max_q, max_u = -1, fill(1.0,n), fill(1,n)
    for q ∈ ProgressBar(LexiIter(probabilities, n))
        for u ∈ Base.product(ntuple(_ -> utilities, n)...)
            # Compute all utilities
            compute_utilities!(single_utils, pair_utils, triple_utils, pools, pairs, triples, q, u)
            # Compute maximum achievable welfare for budgets B ∈ {1,2,3}
            w1 = maximum(single_utils)
            w2 = maximum(pair_utils)
            w3 = maximum(triple_utils)
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
n = 5
@profview faster_search(n, probabilities, utilities)



### The following is a modification to the naive.jl implementation.
### It's currently much slower, for some reason.

markdown"""
TODO:
* Improve the way we iterate over q and u combinations (fully exploit symmetries)
* Extend single computation of utilities / welfares beyond pools to pairs and triples
* Re-introduce multi-threading

# Improve iteration
Suppose $n=3$ with $q=(0.5,0.5,1.0)$. We can assume without loss of generality that the
probabilities are weakly increasing. Then $u=(1,2,3)$ and $u=(2,1,3)$ are essentially
the same (just relabel people 1 and 2). So, for every subgroup of people with the same
probability, we can assume that the utilities are (weakly) increasing. This is
equivalent to using `LexiIter(utilities, k)`, where $k$ is the size of the subgroup.
Note that subgroups are contiguous subsets of $1, \ldots, n$.

This implies a three-step process for iterating over all q and u.
1. Iterate over all ways to sum $n$ as $n = n_1 + \cdots + n_k$ for some $1 \leq k \leq n$.
Each $n_i$ denotes the size of a (contiguous) subgroup.

2. Iterate over `LexiIter(probabilities, k)` to get $(q_1, \ldots, q_k)$, and assign the first
$n_1$ people the probability $q_1$, the next $n_2$ people the probability $q_2$, etc.

3. Iterate over the product of `LexiIter(utilities, n_i)` for each $i \in [k]$.
"""


function improved_search(n::Int, probabilities, utilities)
    pools = collect(powerset(1:n, 1, n)) :: Vector{Vector{Int}}
    allocations = [disjoint_allocations(pools, B) for B ∈ 1:3]
    # Preallocate vectors for computations in loops
    pool_utilities = similar(pools, Float64)  # create an array of Float64s with same length as pools
    w = fill(0.0, 3)  # for storing maximum welfares for budgets 1, 2, 3
    max_diff, max_q, max_u = -1, fill(1.0,n), fill(1,n)
    ends = fill(0, n+1)
    q = fill(0.0, n)
    u = fill(1, n)
    # Start iteration logic for q
    for groupsizes ∈ integer_partitions(n)
        k = length(groupsizes)  # number of subgroups
        ends[2:k+1] .= cumsum(groupsizes)  # index k+1 is where group k ends
        for qs ∈ LexiIter(probabilities, k)
            @inline for i in 1:k; q[ends[i]+1:ends[i+1]] .= qs[i]; end  # set probability for each subgroup
            # Start iteration logic for utilities u
            for us in Base.product(ntuple(i->LexiIter(utilities, groupsizes[i]), k)...)
                @inline for i in 1:k, j in 1:groupsizes[i]; u[ends[i]+j] = us[i][j]; end  # set utilities for all people
                # Compute utility achievable by each pool
                compute_pool_utilities!(pool_utilities, pools, q, u)
                # Compute maximum achievable welfare for budgets B ∈ {1,2,3}
                for B ∈ 1:3
                    w[B] = maximise_welfare(q, u, allocations[B], pool_utilities)
                end
                # Compute marginal return difference, and keep track of solution if it's the current maximum
                diff = w[1] - 2*w[2] + w[3]
                if diff > max_diff + 1e-10
                    max_diff = diff
                    max_q .= q
                    max_u .= u
                end     
            end
        end
    end
    return 
end


probabilities = 0.1:0.1:1.0
utilities = 1:10
n = 4
@profview improved_search(n, probabilities, utilities)


