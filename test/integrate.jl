@testset "integrate" begin
    @testset "integrate_gauss_legendre" begin
        demo_problems = demo_problem.(eachindex(demo_problems_definitions))

        for (f, (a, b), _) in demo_problems
            integral_a_b = ArbExtras.integrate_gauss_legendre(f, Arb(a), Arb(b))
            integral_b_a = ArbExtras.integrate_gauss_legendre(f, Arb(b), Arb(a))
            @test isequal(integral_a_b, -integral_b_a)
            @test Arblib.overlaps(integral_a_b, real(Arblib.integrate(f, a, b)))
        end

        @test Arblib.overlaps(
            ArbExtras.integrate_gauss_legendre(sin, Arb(0), Arb(1)),
            1 - cos(Arb(1)),
        )

        @test Arblib.overlaps(
            ArbExtras.integrate_gauss_legendre(x -> x^3, Arb(0), Arb(1)),
            Arb(1 // 4),
        )

        @test Arblib.overlaps(
            ArbExtras.integrate_gauss_legendre(x -> 1 + x + x^2 + x^3, Arb(1), Arb(3)),
            Arb(104 // 3),
        )

        @test Arblib.overlaps(
            ArbExtras.integrate_gauss_legendre(sin ∘ exp, Arb(0), Arb(4)),
            real(Arblib.integrate(sin ∘ exp, 0, 4)),
        )

        @test isnan(ArbExtras.integrate_gauss_legendre(inv, Arb(-1), Arb(1)))

        # Wide endpoints
        @test Arblib.overlaps(
            ArbExtras.integrate_gauss_legendre(sin, Arb((-1, 0)), Arb((1, 2))),
            cos(Arb((-1, 0))) - cos(Arb((1, 2))),
        )
        # Non-finite endpoints
        @test isnan(ArbExtras.integrate_gauss_legendre(sin, Arb(0), Arb(NaN)))
        @test isnan(ArbExtras.integrate_gauss_legendre(sin, Arb(NaN), Arb(0)))
        @test isnan(ArbExtras.integrate_gauss_legendre(sin, Arb(0), Arb(Inf)))
        @test isnan(ArbExtras.integrate_gauss_legendre(sin, Arb(Inf), Arb(0)))

        # TODO: Test on randomly generated polynomials
    end

    @testset "integrate" begin
        demo_problems = demo_problem.(eachindex(demo_problems_definitions))

        for (f, (a, b), _) in demo_problems
            @test Arblib.overlaps(
                ArbExtras.integrate(f, Arb(a), Arb(b), rtol = Arb("1e-5")),
                real(Arblib.integrate(f, a, b)),
            )
        end

        @test Arblib.overlaps(ArbExtras.integrate(sin, Arb(0), Arb(1)), 1 - cos(Arb(1)))

        @test Arblib.overlaps(ArbExtras.integrate(x -> x^3, Arb(0), Arb(1)), Arb(1 // 4))

        @test Arblib.overlaps(
            ArbExtras.integrate(sin ∘ exp, Arb(0), Arb(4)),
            real(Arblib.integrate(sin ∘ exp, 0, 4)),
        )

        @test isnan(ArbExtras.integrate(inv, Arb(-1), Arb(1)))

        # depth_start > 0
        @test Arblib.overlaps(
            ArbExtras.integrate(sin, Arb(0), Arb(1), depth_start = 4),
            1 - cos(Arb(1)),
        )

        @test Arblib.overlaps(
            ArbExtras.integrate(x -> x^3, Arb(0), Arb(1), depth_start = 5),
            Arb(1 // 4),
        )

        # Non-finite endpoints
        @test isnan(ArbExtras.integrate(sin, Arb(0), Arb(NaN)))
        @test isnan(ArbExtras.integrate(sin, Arb(NaN), Arb(0)))
        @test isnan(ArbExtras.integrate(sin, Arb(0), Arb(Inf)))
        @test isnan(ArbExtras.integrate(sin, Arb(Inf), Arb(0)))
    end
end
