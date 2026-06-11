# Driving a set of solves into the store, with the MOSEK/Gurobi threading split.
# (ProgressMeter is `using`-ed by the module.)

# An algorithm is a NamedTuple `(name, fn, args)`. MOSEK (via greedy's conic
# subroutine) aborts when called from a Julia worker thread on Linux aarch64, so
# greedy solves run on the main thread; Gurobi-based algorithms run in parallel.
main_thread_only(alg) = alg.fn === greedy

"Solve one work item `(pop_hash, pop, T, G, alg)` and record it in the store."
function store_solve!(db, (h, pop, T, G, alg))
    start = Dates.now()
    w, _, guarantee = alg.fn(pop; T=T, G=G, alg.args...)
    record_solve!(db, h, T, G, alg, w, guarantee, Dates.value(Dates.now() - start))
end

"""
    run_experiments(db, algs, populations, budgets, poolsizes; experiment, multithread=false)

Solve every (population, budget, poolsize, algorithm) combination not already in
the store, recording each result. Gurobi solves run in parallel when
`multithread`; MOSEK (greedy) solves always run serially on the main thread.
"""
function run_experiments(db, algs, populations, budgets, poolsizes; experiment, multithread=false)
    hashes = [record_population!(db, pop; experiment, pop_index=i)
              for (i, pop) in enumerate(populations)]
    done = solved_keys(db)
    todo(h, T, G, alg) = (h, T, G, String(alg.name), param_key(alg)) ∉ done

    work = [(h, pop, T, G, alg)
            for (h, pop) in zip(hashes, populations)
            for T in budgets for G in poolsizes for alg in algs if todo(h, T, G, alg)]
    nskip = length(hashes) * length(budgets) * length(poolsizes) * length(algs) - length(work)
    nskip > 0 && println("  ($(nskip) solves already in store, skipping)")

    progress = Progress(length(work); desc="Solving ", showspeed=true)
    label(w) = "$(w[5].name) G=$(w[4]) B=$(w[3]) "
    par = filter(w -> !main_thread_only(w[5]), work)
    serial = filter(w -> main_thread_only(w[5]), work)

    if multithread && !isempty(par)
        # :greedy scheduling pulls cells one at a time, so no thread idles on a slow solve.
        Threads.@threads :greedy for w in par
            store_solve!(db, w); next!(progress; desc=label(w))
        end
    else
        for w in par
            update!(progress; desc=label(w)); store_solve!(db, w); next!(progress; desc=label(w))
        end
    end
    for w in serial  # MOSEK: main thread only
        update!(progress; desc=label(w)); store_solve!(db, w); next!(progress; desc=label(w))
    end
    return nothing
end
