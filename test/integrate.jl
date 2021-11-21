@testset "integrate_gauss_legendre" begin
    # TODO: Test on randomly generated polynomials
    # TODO: Test with alternating endpoints
    # TODO: Test with wide endpoints
    # TODO: Test with non-finite endpoints

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
end

@testset "integrate" begin
    # TODO: Add more tests

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
end
