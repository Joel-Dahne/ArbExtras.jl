@testset "bounded_by" begin
    demo_problems = demo_problem.(1:length(demo_problems_definitions))

    @testset "test problems" begin
        i = 1
        for (f, (a, b), _, fmax, _, fabsmax) in demo_problems
            # We check that the bound is satisfied for a value
            # slightly larger than the maximum but not satisfied for a
            # value slightly smaller. What slightly means depends on
            # the degree we use, sqrt(eps(fmax)) for the default and
            # abs(fmax) / 1000 for degree = -1.
            @test ArbExtras.bounded_by(f, a, b, ubound(fmax + sqrt(eps(fmax))))
            @test !ArbExtras.bounded_by(f, a, b, lbound(fmax - sqrt(eps(fmax))))

            @test ArbExtras.bounded_by(f, a, b, ubound(fmax + abs(fmax) / 100), degree = -1)
            @test !ArbExtras.bounded_by(
                f,
                a,
                b,
                lbound(fmax - abs(fmax) / 100),
                degree = -1,
            )


            @test ArbExtras.bounded_by(
                f,
                a,
                b,
                ubound(fabsmax + sqrt(eps(fabsmax))),
                abs_value = true,
            )
            @test !ArbExtras.bounded_by(
                f,
                a,
                b,
                lbound(fabsmax - sqrt(eps(fabsmax))),
                abs_value = true,
            )

            @test ArbExtras.bounded_by(
                f,
                a,
                b,
                ubound(fabsmax + abs(fabsmax) / 100),
                abs_value = true,
                degree = -1,
            )
            @test !ArbExtras.bounded_by(
                f,
                a,
                b,
                lbound(fabsmax - abs(fabsmax) / 100),
                abs_value = true,
                degree = -1,
            )
        end
    end

    @testset "special cases" begin
        # Thin interval
        @test ArbExtras.bounded_by(x -> 2x, Arf(2), Arf(2), Arf(4))
        @test !ArbExtras.bounded_by(x -> 2x, Arf(2), Arf(2), Arf(3))
        @test ArbExtras.bounded_by(x -> -2x, Arf(2), Arf(2), Arf(4), abs_value = true)
        @test !ArbExtras.bounded_by(x -> -2x, Arf(2), Arf(2), Arf(3), abs_value = true)

        # Non-valid endpoints
        @test_throws ArgumentError ArbExtras.bounded_by(identity, Arf(2), Arf(1), Arf(0))

        # Threading enabled
        @test ArbExtras.bounded_by(x -> 2x, Arf(2), Arf(3), Arf(6.1), threaded = true)
        @test !ArbExtras.bounded_by(x -> 2x, Arf(2), Arf(3), Arf(5.9), threaded = true)

        # verbose = true
        @test_logs (:warn, "reached maximum depth 20") match_mode = :any ArbExtras.bounded_by(
            acos,
            Arf(0),
            Arf(0.5),
            ubound(Arb(π) / 2),
            degree = -1,
            verbose = true,
        )

        @test_logs (:warn, "reached maximum number of evaluations 101 >= 100") match_mode =
            :any ArbExtras.bounded_by(
            sin ∘ exp,
            Arf(0),
            Arf(8),
            Arf(0.999),
            degree = -1,
            maxevals = 100,
            verbose = true,
        )
    end


end
