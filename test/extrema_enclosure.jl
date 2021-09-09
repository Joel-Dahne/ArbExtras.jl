@testset "extrema_enclosure" begin
    @testset "test problems" begin
        demo_problems = demo_problem.(1:length(demo_problems_definitions))

        for (f, (a, b), min, max, absmin, absmax) in demo_problems
            for abs_value in (false, true)
                # Direct evaluation
                min1, max1 = ArbExtras.extrema_enclosure(f, a, b, degree = -1; abs_value)
                min2 = ArbExtras.minimum_enclosure(f, a, b, degree = -1; abs_value)
                max2 = ArbExtras.maximum_enclosure(f, a, b, degree = -1; abs_value)

                # Taylor evaluation
                min3, max3 = ArbExtras.extrema_enclosure(f, a, b; abs_value)
                min4 = ArbExtras.minimum_enclosure(f, a, b; abs_value)
                max4 = ArbExtras.maximum_enclosure(f, a, b; abs_value)

                if abs_value
                    @test Arblib.overlaps(min1, absmin)
                    @test Arblib.overlaps(min2, absmin)
                    @test Arblib.overlaps(min3, absmin)
                    @test Arblib.overlaps(min4, absmin)

                    @test Arblib.overlaps(max1, absmax)
                    @test Arblib.overlaps(max2, absmax)
                    @test Arblib.overlaps(max3, absmax)
                    @test Arblib.overlaps(max4, absmax)
                else
                    @test Arblib.overlaps(min1, min)
                    @test Arblib.overlaps(min2, min)
                    @test Arblib.overlaps(min3, min)
                    @test Arblib.overlaps(min4, min)

                    @test Arblib.overlaps(max1, max)
                    @test Arblib.overlaps(max2, max)
                    @test Arblib.overlaps(max3, max)
                    @test Arblib.overlaps(max4, max)
                end

                # Check that the enclosure satisfies the relative error
                if abs_value && Arblib.contains_zero(min3)
                    # This case is not handled very well at the moment
                    # so it might not return exactly zero for the
                    # minimum, therefore we have to check the absolute
                    # tolerance instead of the relative one.
                    @test Arblib.radius(Arb, min3) <= sqrt(eps(Arb))
                else
                    @test Arblib.radius(Arb, min3) <= sqrt(eps(Arb)) * abs(min3)
                end
                @test Arblib.radius(Arb, max3) <= sqrt(eps(Arb)) * abs(max3)
                if abs_value && Arblib.contains_zero(min4)
                    # This case is not handled very well at the moment
                    # so it might not return exactly zero for the
                    # minimum, therefore we have to check the absolute
                    # tolerance instead of the relative one.
                    @test Arblib.radius(Arb, min4) <= sqrt(eps(Arb))
                else
                    @test Arblib.radius(Arb, min4) <= sqrt(eps(Arb)) * abs(min4)
                end
                @test Arblib.radius(Arb, max4) <= sqrt(eps(Arb)) * abs(max4)
            end
        end
    end

    @testset "special cases" begin
        # Thin interval
        @test ArbExtras.extrema_enclosure(x -> -2x, Arf(2), Arf(2)) == (Arb(-4), Arb(-4))
        @test ArbExtras.minimum_enclosure(x -> -2x, Arf(2), Arf(2)) == Arb(-4)
        @test ArbExtras.maximum_enclosure(x -> -2x, Arf(2), Arf(2)) == Arb(-4)
        @test ArbExtras.extrema_enclosure(x -> -2x, Arf(2), Arf(2), abs_value = true) ==
              (Arb(4), Arb(4))
        @test ArbExtras.minimum_enclosure(x -> -2x, Arf(2), Arf(2), abs_value = true) ==
              Arb(4)
        @test ArbExtras.maximum_enclosure(x -> -2x, Arf(2), Arf(2), abs_value = true) ==
              Arb(4)

        # Non-ordered endpoints
        @test_throws ArgumentError ArbExtras.extrema_enclosure(identity, Arf(2), Arf(1))
        @test_throws ArgumentError ArbExtras.minimum_enclosure(identity, Arf(2), Arf(1))
        @test_throws ArgumentError ArbExtras.maximum_enclosure(identity, Arf(2), Arf(1))

        # Infinite endpoints
        @test_throws ArgumentError ArbExtras.extrema_enclosure(identity, Arf(2), Arf(Inf))
        @test_throws ArgumentError ArbExtras.minimum_enclosure(identity, Arf(2), Arf(Inf))
        @test_throws ArgumentError ArbExtras.maximum_enclosure(identity, Arf(2), Arf(Inf))

        # Threading enabled
        @test Arblib.isequal(
            ArbExtras.extrema_enclosure(
                x -> -2x,
                Arf(2),
                Arf(3),
                degree = -1,
                threaded = true,
            ),
            ArbExtras.extrema_enclosure(x -> -2x, Arf(2), Arf(3), degree = -1),
        )
        @test Arblib.isequal(
            ArbExtras.extrema_enclosure(x -> -2x, Arf(2), Arf(3), threaded = true),
            ArbExtras.extrema_enclosure(x -> -2x, Arf(2), Arf(3)),
        )
        @test Arblib.isequal(
            ArbExtras.minimum_enclosure(
                x -> -2x,
                Arf(2),
                Arf(3),
                degree = -1,
                threaded = true,
            ),
            ArbExtras.minimum_enclosure(x -> -2x, Arf(2), Arf(3), degree = -1),
        )
        @test Arblib.isequal(
            ArbExtras.minimum_enclosure(x -> -2x, Arf(2), Arf(3), threaded = true),
            ArbExtras.minimum_enclosure(x -> -2x, Arf(2), Arf(3)),
        )
        @test Arblib.isequal(
            ArbExtras.maximum_enclosure(
                x -> -2x,
                Arf(2),
                Arf(3),
                degree = -1,
                threaded = true,
            ),
            ArbExtras.maximum_enclosure(x -> -2x, Arf(2), Arf(3), degree = -1),
        )
        @test Arblib.isequal(
            ArbExtras.maximum_enclosure(x -> -2x, Arf(2), Arf(3), threaded = true),
            ArbExtras.maximum_enclosure(x -> -2x, Arf(2), Arf(3)),
        )

        # Verbose = true
        @test_logs (:info, "reached maximum depth 20") match_mode = :any ArbExtras.extrema_enclosure(
            inv,
            Arf(0),
            Arf(1),
            verbose = true,
        )
        @test_logs (:info, "reached maximum depth 20") match_mode = :any ArbExtras.minimum_enclosure(
            inv,
            Arf(0),
            Arf(1),
            verbose = true,
        )
        @test_logs (:info, "reached maximum depth 20") match_mode = :any ArbExtras.maximum_enclosure(
            inv,
            Arf(0),
            Arf(1),
            verbose = true,
        )

        @test_logs (:info, "reached maximum number of evaluations 117 >= 100") match_mode =
            :any ArbExtras.extrema_enclosure(
            sin ∘ exp,
            Arf(0),
            Arf(10),
            maxevals = 100,
            verbose = true,
        )
        @test_logs (:info, "reached maximum number of evaluations 107 >= 100") match_mode =
            :any ArbExtras.minimum_enclosure(
            sin ∘ exp,
            Arf(0),
            Arf(10),
            maxevals = 100,
            verbose = true,
        )
        @test_logs (:info, "reached maximum number of evaluations 111 >= 100") match_mode =
            :any ArbExtras.maximum_enclosure(
            sin ∘ exp,
            Arf(0),
            Arf(10),
            maxevals = 100,
            verbose = true,
        )
    end
end
