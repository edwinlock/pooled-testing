@testset "datastore" begin

    alg(name; args...) = (name=name, fn=identity, args=Dict(args...))

    @testset "population_hash is content-based and order-independent" begin
        a = Dict(1 => (0.9, 10), 2 => (0.8, 5))
        b = Dict(2 => (0.8, 5), 1 => (0.9, 10))   # same content, different order
        c = Dict(1 => (0.9, 10), 2 => (0.8, 6))   # one utility differs
        @test population_hash(a) == population_hash(b)
        @test population_hash(a) != population_hash(c)
        @test length(population_hash(a)) == 16
    end

    @testset "param_key captures solve-affecting args, ignores verbose" begin
        @test param_key(alg(:approx; K=15, verbose=false)) ==
              param_key(alg(:approx; K=15, verbose=true))   # verbose excluded
        @test param_key(alg(:approx; K=15)) != param_key(alg(:approx; K=20))
        @test occursin("mipgap=$(GUROBI_MIPGAP)", param_key(alg(:greedy)))
    end

    @testset "record + retrieve + idempotency" begin
        db = open_store(mktempdir())
        pop = Dict(1 => (0.9, 10), 2 => (0.8, 5))
        h = record_population!(db, pop; experiment=3, pop_index=1)
        @test h == population_hash(pop)

        a = alg(:approx; K=15, verbose=false)
        record_solve!(db, h, 6, 5, a, 100.0, 0.5, 1234)
        record_solve!(db, h, 6, 5, a, 999.0, 0.5, 1234)   # same key ⇒ replace, not duplicate

        keys = solved_keys(db)
        @test (h, 6, 5, "approx", param_key(a)) in keys
        @test length(keys) == 1                           # idempotent: one row, not two
    end

    @testset "solves are reused across experiments (same content)" begin
        db = open_store(mktempdir())
        pop = Dict(1 => (0.9, 10), 2 => (0.8, 5))
        h1 = record_population!(db, pop; experiment=3, pop_index=1)
        h2 = record_population!(db, pop; experiment=4, pop_index=7)  # same content
        @test h1 == h2                                    # ⇒ same solve key
        a = alg(:approx; K=15)
        record_solve!(db, h1, 6, 5, a, 100.0, 0.5, 10)
        @test (h2, 6, 5, "approx", param_key(a)) in solved_keys(db)  # exp 4 reuses exp 3's solve
    end

    @testset "load_solves joins population metadata" begin
        rootdir = mktempdir()
        db = open_store(rootdir)
        h = record_population!(db, Dict(j => (0.8, 5) for j in 1:9); experiment=3, pop_index=2)
        record_solve!(db, h, 4, 5, alg(:greedy), 50.0, 0.0, 20)
        df = load_solves(db)   # reuse the open handle (avoids a second connection)
        @test nrow(df) == 1
        @test df.experiment[1] == 3 && df.pop_index[1] == 2 && df.n[1] == 9
        @test "instance_type" in names(df) && "gurobi_version" in names(df)
    end

end
