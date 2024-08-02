@testset "SpecialFunctions" begin
    @testset "besselj and bessely" begin
        # Check recurrence relations
        # besselj(ν, z) = z / 2ν * (besselj(ν - 1, z) + besselj(ν + 1, z))
        # and same for bessely
        rng = MersenneTwister(42)
        setprecision(Arb, 64) do
            for _ = 1:1000
                deg = rand(rng, 0:10)
                ν = Arb(20rand(rng) - 10)
                z = ArbSeries(20rand(rng, deg + 1) .- 10)

                lhsj = besselj(ν, z)
                rhsj = z / 2ν * (besselj(ν - 1, z) + besselj(ν + 1, z))
                @test Arblib.degree(lhsj) == Arblib.degree(rhsj) == deg
                @test Arblib.overlaps.(lhsj, rhsj)

                lhsy = bessely(ν, z)
                rhsy = z / 2ν * (bessely(ν - 1, z) + bessely(ν + 1, z))
                @test Arblib.degree(lhsy) == Arblib.degree(rhsy) == deg
                @test Arblib.overlaps(lhsy, rhsy)
            end
        end

        z = ArbSeries([1, 2, 3, 4, 5])
        @test isequal(besselj0(z), besselj(zero(Arb), z))
        @test isequal(besselj1(z), besselj(one(Arb), z))
        @test isequal(bessely0(z), bessely(zero(Arb), z))
        @test isequal(bessely1(z), bessely(one(Arb), z))
    end
end
