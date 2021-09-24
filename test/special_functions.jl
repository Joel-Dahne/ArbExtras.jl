@testset "SpecialFunctions" begin
    @testset "besselj" begin
        # Check recurrence relations
        # besselj(ν, z) = z / 2ν * (besselj(ν - 1, z) + besselj(ν + 1, z))
        rng = MersenneTwister(42)
        setprecision(Arb, 64) do
            for _ = 1:1000
                deg = rand(rng, 0:10)
                ν = Arb(20rand(rng) - 10)
                z = ArbSeries(20rand(rng, deg + 1) .- 10)

                lhs = besselj(ν, z)
                rhs = z / 2ν * (besselj(ν - 1, z) + besselj(ν + 1, z))
                @test Arblib.degree(lhs) == Arblib.degree(rhs) == deg
                @test all(Arblib.overlaps.(lhs[:], rhs[:]))
            end
        end

        z = ArbSeries([1, 2, 3, 4, 5])
        @test isequal(besselj0(z), besselj(zero(Arb), z))
        @test isequal(besselj1(z), besselj(one(Arb), z))
    end
end
