@testset "temporary" begin
    @testset "iscpx" begin
        @test !ArbExtras.iscpx(ArbPoly())
        @test !ArbExtras.iscpx(ArbSeries())
        @test !ArbExtras.iscpx(AcbPoly())
        @test !ArbExtras.iscpx(AcbSeries())

        @test !ArbExtras.iscpx(ArbPoly(5))
        @test !ArbExtras.iscpx(ArbSeries(5))
        @test !ArbExtras.iscpx(AcbPoly(5))
        @test !ArbExtras.iscpx(AcbSeries(5))

        @test ArbExtras.iscpx(ArbPoly((5, 1)))
        @test ArbExtras.iscpx(ArbSeries((5, 1)))
        @test ArbExtras.iscpx(AcbPoly((5, 1)))
        @test ArbExtras.iscpx(AcbSeries((5, 1)))

        @test !ArbExtras.iscpx(ArbPoly((5, 1, 1)))
        @test !ArbExtras.iscpx(ArbSeries((5, 1, 1)))
        @test !ArbExtras.iscpx(AcbPoly((5, 1, 1)))
        @test !ArbExtras.iscpx(AcbSeries((5, 1, 1)))

        @test !ArbExtras.iscpx(ArbPoly((5, 2)))
        @test !ArbExtras.iscpx(ArbSeries((5, 2)))
        @test !ArbExtras.iscpx(AcbPoly((5, 2)))
        @test !ArbExtras.iscpx(AcbSeries((5, 2)))

        @test ArbExtras.iscpx(ArbSeries((5, 1, 0)))
        @test ArbExtras.iscpx(AcbPoly((5, 1, 0)))

        @test !ArbExtras.iscpx(ArbSeries((5, 0)))
        @test !ArbExtras.iscpx(AcbPoly((5, 0)))
    end

    @testset "compose_zero" begin
        for coeffs_p in ([1, 2, 3], [5, 2], [])
            for coeffs_q in ([5, 1], [5, 2], [1, 2, 3, 5, 7], [0, 2, 3, 5, 7], [])

                coeffs0_q = copy(coeffs_q)
                if !isempty(coeffs_q)
                    coeffs0_q[1] = 0
                end

                @test Arblib.overlaps(
                    Arblib.compose(ArbPoly(coeffs_p), ArbPoly(coeffs0_q)),
                    ArbExtras.compose_zero(ArbPoly(coeffs_p), ArbPoly(coeffs_q)),
                )
            end
        end
    end
end
