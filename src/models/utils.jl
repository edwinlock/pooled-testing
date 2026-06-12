"""
Generate a random instance of the test allocation problem.
Returns two vectors of length n whose ith entries denote the
probability of being healthy and utility for the i-th
person.

If `utils` [`probs`] is a vector or a range, entries are chosen uniformly at
random. Otherwise, it is a probability distribution, and entries are sampled
appropriately.

Utilities are rounded to the nearest integer.
"""
function generate_instance(n, probs, utils)
	# generate random utilities and probabilities of being healthy
    pop = Population{Int}()
    for i in 1:n
        pop[i] = (rand(probs), round(rand(utils)))
    end
	return pop
end


"""
Scale and round utilities of population `pop` so that they are integers between 1 and `upper`.
"""
function scale_utilities(pop::Population, upper)
    isempty(pop) && return pop
    max_util = maximum(val[2] for (_, val) in pop)
    factor = max_util <= upper ? 1 : upper / max_util  #  no need to scale if max_util <= upper
    output = empty(pop)
    for (key, val) in pop
        output[key] = (val[1], clamp(round(Int, factor*val[2]), 1, upper))
    end
    return output
end


"""
Remove people with zero (or negative) utilities from the population.
Returns new population dict.
"""
remove_zeros!(pop::Population) = filter!(x->x.second[1]>0&&x.second[2] > 0, pop)


function pop2vec(pop)
    keylist = collect(keys(pop))
    q = [pop[k][1] for k in keylist]
    u = [pop[k][2] for k in keylist]
    return q, u, keylist
end



"""
Cluster input population `pop' by grouping together everyone with the same
(q, u) tuple. Keys in `clusters` are 
"""
function pop2clusters(pop::Population)
    clusters = DefaultDict(Vector)
    for (key, val) in pop
        push!(clusters[val], key)
    end
    return clusters
end


"""
Convert a cluster dict to vectors.
"""
function cluster2vec(d)
    keylist = collect(keys(d))
    q = [k[1] for k in keylist]
    u = [k[2] for k in keylist]
    n = [length(d[k]) for k in keylist]
    return q, u, n, keylist
end


"""
Given `pools` as lists of entries that are indices to the `clusters`,
replace indices with representatives from the correct clusters using `keys`.
"""
function uncluster(pools, clusters, keylist)
    new_pools = []
    for pool in pools
        new_pool = []
        for i in pool
            person = pop!(clusters[keylist[i]])
            push!(new_pool, person)
        end
        push!(new_pools, new_pool)
    end
    return new_pools
end


"""Retrieve welfare and pools from solution of model with variables x."""
function retrieve(m, x, T, C)
    welfare = objective_value(m)
    v = Int.(round.(value.(x)))
    pools = [vcat([fill(i, v[t,i]) for i in 1:C]...) for t in 1:T]
    return welfare, pools
end


function welfare(pools, pop)
    w = 0
    for s in powerset(pools, 1)
        utilsum = sum(pop[k][2] for k in intersect(s...); init=0)
        healthprod = prod(pop[k][1] for k in union(s...); init=1)
        w += (-1)^(length(s)+1)*utilsum*healthprod
    end
    return w
end


"""
Set start values for x variables in model m based on `pools`.
"""
function start_vals(m, pools, keylist)
    x = m[:x]
    set_start_value.(x, 0)  # exclude everyone
    for (i, pool) in enumerate(pools)
        i > size(x, 1) && break
        for person in pool
            j = findfirst(==(person), keylist)  # find index of person with id `person`
            isnothing(j) && error("Person $(person) is not in keylist.")
            set_start_value(x[i,j], start_value(x[i,j]) + 1)
        end
    end
end

"Index of the PWL segment (per breakpoints `c`) containing `l`."
function exp_segment(l, c)
    for k in 1:(length(c)-1)
        if c[k] - 1e-9 <= l <= c[k+1] + 1e-9
            return k
        end
    end
    return l < c[1] ? 1 : length(c)-1
end

"""
Set a fuller clustered MIP start from pools of person ids.

In addition to `x`, this fills the deterministic variables implied by each
started pool: `z`, `zind`, `y`, `l`, `lind`, `v`, and a feasible `w`.
"""
function full_clustered_start_vals(m, pools, pop, keylist, q, u; G=5, K=15)
    x = m[:x]
    z = m[:z]
    y = m[:y]
    l = m[:l]
    w = m[:w]
    zind = m[:zind]
    lind = m[:lind]
    v = m[:v]

    A, B = exp_domain(q, u, G)
    a, b, c, _ = upper_bounds(A, B, K)

    for (t, pool) in enumerate(pools)
        t > size(x, 1) && break
        counts = zeros(Int, length(keylist))
        log_health = 0.0
        utility_sum = 0

        for person in pool
            cluster_key = pop[person]
            j = findfirst(==(cluster_key), keylist)
            isnothing(j) && error("Cluster $(cluster_key) for person $(person) is not in keylist.")
            counts[j] += 1
            log_health += log(cluster_key[1])
            utility_sum += cluster_key[2]
        end

        for j in eachindex(keylist)
            set_start_value(x[t,j], counts[j])
        end

        yval = log(utility_sum)
        lval = yval + log_health
        active_segment = exp_segment(lval, c)
        set_start_value(z[t], utility_sum)
        set_start_value(y[t], yval)
        set_start_value(l[t], lval)
        set_start_value(w[t], a[active_segment] * lval + b[active_segment])

        for k in axes(zind, 2)
            set_start_value(zind[t,k], k == utility_sum ? 1 : 0)
        end

        for k in 1:K
            set_start_value(lind[t,k], k == active_segment ? 1 : 0)
            set_start_value(v[t,k], k == active_segment ? lval : 0)
        end
    end
end
