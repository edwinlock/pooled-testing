"""
Generate a random instance of the test allocation problem.
Returns two vectors of length n whose ith entries denote the
probability of being healthy and utility for the i-th
person.

If `utils` [`probs`] is a vector, entries are chosen uniformly at
random. Otherwise, it is a probability distribution.

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
        output[key] = (val[1], Int(round(factor*val[2])))
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
        for person in pool
            j = findall(x->x==2, keylist)[1]  # find index of person with id `person`
            set_start_value(x[i,j], start_value(x[i,j])+1)
        end
    end
end