@testset "fx_div_x" begin
    fs = [sin, tan, x -> x^2, x -> exp(x) - 1, x -> sin(x) - x, x -> sin(x + π)]
    orders = [1, 1, 1, 1, 2, 1]
    forces = [false, false, false, false, false, true]

    @testset "x::Arb" begin
        for (f, order, force) in zip(fs, orders, forces)
            for radius in Arb(2) .^ (-10:1)
                x = Arblib.add_error!(zero(Arb), radius)
                if isfinite(f(ArbSeries((x, 1)))) # No point to test otherwise
                    enclosures = [
                        ArbExtras.fx_div_x(
                            f,
                            x,
                            order;
                            extra_degree,
                            enclosure_degree,
                            force,
                        ) for extra_degree = 0:2, enclosure_degree = -1:2
                    ]
                    @test all(isfinite, enclosures)
                    for y in range(getinterval(Arb, x)..., 10)
                        exact = f(y) / y^order
                        @test all(Arblib.overlaps.(enclosures, Ref(exact)))
                    end
                end
            end
        end
    end

    @testset "x::Acb" begin
        for (f, order, force) in zip(fs, orders, forces)
            for radius in Arb(2) .^ (-10:1)
                x = Arblib.add_error!(zero(Acb), radius)
                if isfinite(f(AcbSeries((x, 1)))) # No point to test otherwise
                    enclosures = [
                        ArbExtras.fx_div_x(
                            f,
                            x,
                            order;
                            extra_degree,
                            enclosure_degree,
                            force,
                        ) for extra_degree = 0:2, enclosure_degree = -1:2
                    ]
                    @test all(isfinite, enclosures)
                    y_reals = range(getinterval(Arb, real(x))..., 10)
                    y_imags = range(getinterval(Arb, imag(x))..., 10)
                    for y_real in y_reals
                        y1 = Acb(y_real, y_imags[1])
                        y2 = Acb(y_real, y_imags[end])
                        @test all(Arblib.overlaps.(enclosures, Ref(f(y1) / y1^order)))
                        @test all(Arblib.overlaps.(enclosures, Ref(f(y2) / y2^order)))
                    end
                    for y_imag in y_imags
                        y1 = Acb(y_reals[1], y_imag)
                        y2 = Acb(y_reals[end], y_imag)
                        @test all(Arblib.overlaps.(enclosures, Ref(f(y1) / y1^order)))
                        @test all(Arblib.overlaps.(enclosures, Ref(f(y2) / y2^order)))
                    end
                end
            end
        end
    end

    @testset "x::ArbSeries" begin
        # Pure derivative
        for (f, order, force) in zip(fs, orders, forces)
            for radius in Arb(2) .^ (-10:1)
                for degree = 0:4
                    x = ArbSeries((Arblib.add_error!(zero(Arb), radius), 1); degree)
                    if isfinite(f(ArbSeries((x[0], 1)))) # No point to test otherwise
                        enclosures = [
                            ArbExtras.fx_div_x(
                                f,
                                x,
                                order;
                                extra_degree,
                                enclosure_degree,
                                force,
                            ) for extra_degree = 0:2, enclosure_degree = -1:2
                        ]
                        @test all(isfinite, enclosures)
                        for y in range(getinterval(Arb, x[0])..., 10)
                            y_s = ArbSeries(x)
                            y_s[0] = y
                            exact = f(y_s) / y_s^order
                            @test all(Arblib.overlaps.(enclosures, Ref(exact)))
                        end
                    end
                end
            end
        end

        # General series
        for (f, order, force) in zip(fs, orders, forces)
            for radius in Arb(2) .^ (-10:1)
                for degree = 0:4
                    x = ArbSeries((Arblib.add_error!(zero(Arb), radius), 2, 3); degree)
                    if isfinite(f(ArbSeries((x[0], 1)))) # No point to test otherwise
                        enclosures = [
                            ArbExtras.fx_div_x(
                                f,
                                x,
                                order;
                                extra_degree,
                                enclosure_degree,
                                force,
                            ) for extra_degree = 0:2, enclosure_degree = -1:2
                        ]
                        @test all(isfinite, enclosures)
                        for y in range(getinterval(Arb, x[0])..., 10)
                            y_s = ArbSeries(x)
                            y_s[0] = y
                            exact = f(y_s) / y_s^order
                            @test all(Arblib.overlaps.(enclosures, Ref(exact)))
                        end
                    end
                end
            end
        end
    end

    @testset "x::AcbSeries" begin
        # Pure derivative
        for (f, order, force) in zip(fs, orders, forces)
            for radius in Arb(2) .^ (-10:1)
                for degree = 0:4
                    x = AcbSeries((Arblib.add_error!(zero(Acb), radius), 1); degree)
                    if isfinite(f(AcbSeries((x[0], 1)))) # No point to test otherwise
                        enclosures = [
                            ArbExtras.fx_div_x(
                                f,
                                x,
                                order;
                                extra_degree,
                                enclosure_degree,
                                force,
                            ) for extra_degree = 0:2, enclosure_degree = -1:2
                        ]
                        @test all(isfinite, enclosures)
                        y_reals = range(getinterval(Arb, real(x[0]))..., 10)
                        y_imags = range(getinterval(Arb, imag(x[0]))..., 10)
                        for y_real in y_reals
                            y1 = AcbSeries(x);
                            y1[0] = Acb(y_real, y_imags[1])
                            y2 = AcbSeries(x);
                            y2[0] = Acb(y_real, y_imags[end])
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y1) / y1^order)))
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y2) / y2^order)))
                        end
                        for y_imag in y_imags
                            y1 = AcbSeries(x);
                            y1[0] = Acb(y_reals[1], y_imag)
                            y2 = AcbSeries(x);
                            y2[0] = Acb(y_reals[end], y_imag)
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y1) / y1^order)))
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y2) / y2^order)))
                        end
                    end
                end
            end
        end

        # General series
        for (f, order, force) in zip(fs, orders, forces)
            for radius in Arb(2) .^ (-10:1)
                for degree = 0:4
                    x = AcbSeries((Arblib.add_error!(zero(Acb), radius), 2, 3); degree)
                    if isfinite(f(AcbSeries((x[0], 1)))) # No point to test otherwise
                        enclosures = [
                            ArbExtras.fx_div_x(
                                f,
                                x,
                                order;
                                extra_degree,
                                enclosure_degree,
                                force,
                            ) for extra_degree = 0:2, enclosure_degree = -1:2
                        ]
                        @test all(isfinite, enclosures)
                        y_reals = range(getinterval(Arb, real(x[0]))..., 10)
                        y_imags = range(getinterval(Arb, imag(x[0]))..., 10)
                        for y_real in y_reals
                            y1 = AcbSeries(x);
                            y1[0] = Acb(y_real, y_imags[1])
                            y2 = AcbSeries(x);
                            y2[0] = Acb(y_real, y_imags[end])
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y1) / y1^order)))
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y2) / y2^order)))
                        end
                        for y_imag in y_imags
                            y1 = AcbSeries(x);
                            y1[0] = Acb(y_reals[1], y_imag)
                            y2 = AcbSeries(x);
                            y2[0] = Acb(y_reals[end], y_imag)
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y1) / y1^order)))
                            @test all(Arblib.overlaps.(enclosures, Ref(f(y2) / y2^order)))
                        end
                    end
                end
            end
        end
    end

    @testset "Non-finite" begin
        x_Arb = Arblib.add_error!(zero(Arb), Arb(1))
        x_Acb = Arblib.add_error!(zero(Acb), Arb(1))
        x_ArbSeries = ArbSeries((Arblib.add_error!(zero(Arb), Arb(1)), 1))
        x_AcbSeries = AcbSeries((Arblib.add_error!(zero(Acb), Arb(1)), 1))

        @test !isfinite(ArbExtras.fx_div_x(x -> 1 / x, x_Arb))
        @test !isfinite(ArbExtras.fx_div_x(x -> 1 / x, x_Acb))
        @test !isfinite(ArbExtras.fx_div_x(x -> 1 / x, x_ArbSeries))
        @test !isfinite(ArbExtras.fx_div_x(x -> 1 / x, x_AcbSeries))
    end

    @testset "Zero" begin
        # sin(x) / x → 1 as x → 0
        @test Arblib.overlaps(ArbExtras.fx_div_x(sin, zero(Arb)), Arb(1))
        @test Arblib.overlaps(ArbExtras.fx_div_x(sin, zero(Acb)), Acb(1))
        # sin(x) / x at x = 0 + 1*t gives the series [1, 0]: constant term 1, linear term 0
        @test Arblib.overlaps(
            ArbExtras.fx_div_x(sin, ArbSeries((zero(Arb), 1))),
            ArbSeries((1, 0)),
        )
        @test Arblib.overlaps(
            ArbExtras.fx_div_x(sin, AcbSeries((zero(Acb), 1))),
            AcbSeries((1, 0)),
        )
    end

    @testset "Throws" begin
        x_Arb = Arblib.add_error!(zero(Arb), Arb(1))
        x_Acb = Arblib.add_error!(zero(Acb), Arb(1))
        x_ArbSeries = ArbSeries((Arblib.add_error!(zero(Arb), Arb(1)), 1))
        x_AcbSeries = AcbSeries((Arblib.add_error!(zero(Acb), Arb(1)), 1))

        # x does not contain zero
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, Arb(1))
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, Acb(1))
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, ArbSeries((Arb(1), 1)))
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, AcbSeries((Acb(1), 1)))

        # order < 1
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_Arb, 0)
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_Acb, 0)
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_ArbSeries, 0)
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_AcbSeries, 0)

        # extra_degree < 0
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_Arb; extra_degree = -1)
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_Acb; extra_degree = -1)
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_ArbSeries; extra_degree = -1)
        @test_throws ArgumentError ArbExtras.fx_div_x(sin, x_AcbSeries; extra_degree = -1)

        # force = true but f does not have a zero at 0
        @test_throws ErrorException ArbExtras.fx_div_x(cos, x_Arb; force = true)
        @test_throws ErrorException ArbExtras.fx_div_x(cos, x_Acb; force = true)
        @test_throws ErrorException ArbExtras.fx_div_x(cos, x_ArbSeries; force = true)
        @test_throws ErrorException ArbExtras.fx_div_x(cos, x_AcbSeries; force = true)
    end
end
