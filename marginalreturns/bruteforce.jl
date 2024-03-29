using Combinatorics, ProgressBars

include("iterator.jl")

"""Compute all disjoint allocations for population size n and B tests."""
all_allocations(n, B) = [c for partitions in [collect(partitions(s, B)) for subset in [collect(combinations(1:n, k)) for k in B:n] for s in subset] for c in partitions]

"""Maximal welfare over set of allocations given."""
function opt_welfare(q, u, allocations)
    max_welfare, welfare = -1.0, 0.0
    for a in allocations
        welfare = 0.0  # reset for each allocation
        for pool in a
            welfare += prod(q[i] for i in pool) * sum(u[i] for i in pool)
        end
        max_welfare = max(max_welfare, welfare)
    end
    return max_welfare
end


function search(n::Int, probabilities, utilities)
    allocations = [all_allocations(n, B) for B in 1:3]  # compute all possible allocations
    w = [0.0, 0.0, 0.0]  # pre-allocate vector for storing welfares for budgets 1, 2, 3
    max_diff, max_q, max_u = -1, fill(1.0,n), fill(1,n)
    for q ∈ ProgressBar(LexiIter(probabilities,n))
        for u ∈ Base.product(ntuple(_ -> utilities, n)...)
            for B in 1:3; w[B] = opt_welfare(q, u, allocations[B]); end
            diff = w[1] - 2*w[2] + w[3]
            if diff > max_diff + 1e-10
                max_diff = diff
                max_q .= q
                max_u .= u
            end
        end
    end
    return max_diff, (max_q, max_u)
end

probabilities = 0.1:0.1:1.0
utilities = 1:10
n = 6
search(n, probabilities, utilities)


# function subroutine(n::Int, first_probabilities::Array{Float64}, other_probabilities::Array{Float64}, utilities::Array{Int})
#     allocations = [all_allocations(n, B) for B in 1:3]  # compute all possible allocations
#     w = [0.0, 0.0, 0.0]  # pre-allocate vector for storing welfares for budgets 1, 2, 3
#     max_diff, max_q, max_u = -1, fill(1.0,n), fill(1,n)
#     all_probs = [first_probabilities, ntuple(i->other_probabilities, n-1)...]
#     for q ∈ Base.product(all_probs...)
#         for u ∈ Base.product(ntuple(i->utilities, n)...)
#             for B in 1:3; w[B] = opt_welfare(q, u, allocations[B]); end
#             diff = w[1] - 2*w[2] + w[3]
#             if diff > max_diff + 1e-10
#                 max_diff = diff
#                 max_q .= q
#                 max_u .= u
#             end
#         end
#     end
#     return max_diff, (max_q, max_u)
# end


# function multithreaded(n::Int, probabilities::Array{Float64}, utilities::Array{Int})
#     nthreads = Threads.nthreads()
#     max_diffs = fill(0.0, nthreads)
#     max_sols = fill(( fill(0.0, n), fill(1, n) ), nthreads)
#     part_sizes = ceil(Int, length(probabilities) / nthreads)
#     partition = collect(Iterators.partition(probabilities,part_sizes))
#     @Threads.threads for i in 1:nthreads
#         first_probs = collect(partition[i])
#         max_diffs[i], max_sols[i] = subroutine(n, first_probs, probabilities, utilities)
#     end
#     i = argmax(max_diffs)  # find thread that maximised diff
#     return max_diffs[i], max_sols[i]
# end


# @time multithreaded(n, probabilities, utilities)
