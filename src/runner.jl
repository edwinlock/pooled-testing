# Driving a set of solves into the store. (ProgressMeter is `using`-ed by the
# module.) An algorithm is a NamedTuple `(name, fn, args)`.
#
# NB: MOSEK aborts when called from a Julia worker thread on Linux aarch64. The
# experiments run on x86, where this doesn't apply; if that changes, greedy
# solves must be quarantined to the main thread again (see git history).

"Solve one work item `(experiment, pop_index, pop_hash, pop, T, G, alg)` and record it in the store."
function store_solve!(db, (experiment, pop_index, h, pop, T, G, alg))
    start = Dates.now()
    w, _, guarantee, guarantee_post = alg.fn(pop; T=T, G=G, alg.args...)
    record_solve!(db, experiment, pop_index, h, T, G, alg, w, guarantee, guarantee_post, Dates.value(Dates.now() - start))
end

"""
    run_experiments(db, algs, populations, budgets, poolsizes; experiment, multithread=false)

Solve every (population, budget, poolsize, algorithm) combination not already in
the store, recording each result. Solves run in parallel when `multithread`.
"""
function run_experiments(db, algs, populations, budgets, poolsizes; experiment, multithread=false, force=false)
    hashes = [record_population!(db, pop) for pop in populations]
    done = solved_keys(db)
    # `force` re-solves cells already in the store; record_solve! is INSERT OR
    # REPLACE, so the existing row is overwritten in place (e.g. greedy re-timing).
    todo(i, h, T, G, alg) = force || (experiment, i, h, T, G, String(alg.name), param_key(alg)) ∉ done

    work = [(experiment, i, h, pop, T, G, alg)
            for (i, (h, pop)) in enumerate(zip(hashes, populations))
            for T in budgets for G in poolsizes for alg in algs if todo(i, h, T, G, alg)]
    nskip = length(hashes) * length(budgets) * length(poolsizes) * length(algs) - length(work)
    nskip > 0 && println("  ($(nskip) solves already in store, skipping)")

    progress = Progress(length(work); dt=0.0, desc="Solving ", showspeed=true)
    label(w) = "$(w[7].name) G=$(w[6]) B=$(w[5]) "

    if multithread
        # :greedy scheduling pulls cells one at a time, so no thread idles on a slow solve.
        Threads.@threads :greedy for w in work
            ProgressMeter.update!(progress; desc=label(w))
            store_solve!(db, w)
            ProgressMeter.next!(progress; desc=label(w))
        end
    else
        for w in work
            ProgressMeter.update!(progress; desc=label(w))
            store_solve!(db, w)
            ProgressMeter.next!(progress; desc=label(w))
        end
    end
    return nothing
end
