@testset "analyse: reshaping the store" begin
    include(joinpath(REPO, "analyse.jl"))   # brings in experiment_frame, summary_table

    # Build a small store: 2 populations × 2 budgets × {approx, greedy}.
    rootdir = mktempdir()
    db = open_store(rootdir)
    for pi in 1:2, T in (2, 4)
        # Distinct content per pop_index ⇒ distinct hash (else the store dedupes them).
        h = record_population!(db, Dict(j => (0.8, 5 + pi) for j in 1:5); experiment=3, pop_index=pi)
        record_solve!(db, h, T, 5, (name=:approx, fn=identity, args=Dict(:K=>15)), 100.0+T+pi, 0.5, 1000)
        record_solve!(db, h, T, 5, (name=:greedy, fn=identity, args=Dict()),       100.0+T+pi-0.3, 0.0, 50)
    end

    @testset "experiment_frame pivots to wide" begin
        w = experiment_frame(rootdir, 3; algorder=["approx", "greedy"])
        @test nrow(w) == 4                                    # 2 pops × 2 budgets
        for c in ["approx_welfare","greedy_welfare","approx_time","greedy_time","diff","ratio"]
            @test c in names(w)
        end
        @test all(w.diff .≈ 0.3)                              # approx − greedy
        @test all(w.ratio .> 1)                               # approx ≥ greedy
        @test issorted(w.budget)
    end

    @testset "experiment_frame empty for an unrun experiment" begin
        @test isempty(experiment_frame(rootdir, 99; algorder=["approx","greedy"]))
    end

    @testset "summary_table aggregates per budget" begin
        w = experiment_frame(rootdir, 3; algorder=["approx", "greedy"])
        s = summary_table(w, ["approx", "greedy"])
        @test nrow(s) == 2                                    # one row per budget
        @test "approx welfare" in names(s)
        @test "approx guarantee" in names(s)                  # guarantee surfaced
        # Mean welfare at budget 2 over pops {1,2}: (103 + 104)/2 = 103.5
        row = s[s.budget .== 2, :]
        @test row[1, "approx welfare"] ≈ 103.5
    end
end
