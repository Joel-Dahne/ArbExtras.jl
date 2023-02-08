function SpecialFunctions.besselj(ν::Arblib.ArbOrRef, z::Arblib.ArbSeries)
    deg = Arblib.degree(z)
    res = zero(z)

    x = z[0]
    tmp = zero(x)
    ai = zero(x)

    # ai = besselj(ν, x)
    Arblib.hypgeom_bessel_j!(ai, ν, x)

    res[0] = ai

    if deg >= 1
        # ai = (besselj(ν - 1, x) - besselj(ν + 1, x)) / 2

        # ai = besselj(ν - 1, x)
        Arblib.sub!(tmp, ν, 1)
        Arblib.hypgeom_bessel_j!(ai, tmp, x)

        # tmp = besselj(ν + 1, x)
        # ai += tmp
        Arblib.add!(tmp, ν, 1)
        Arblib.hypgeom_bessel_j!(tmp, tmp, x)
        Arblib.sub!(ai, ai, tmp)

        # ai /= 2
        Arblib.mul_2exp!(ai, ai, -1)

        res[1] = ai

        if deg >= 2
            # ai = (besselj(ν - 2, x) + besselj(ν + 2, x) - 2a0) / 8

            # ai = besselj(ν - 2, x)
            Arblib.sub!(tmp, ν, 2)
            Arblib.hypgeom_bessel_j!(ai, tmp, x)

            # tmp = besselj(ν + 2, x)
            # ai += tmp
            Arblib.add!(tmp, ν, 2)
            Arblib.hypgeom_bessel_j!(tmp, tmp, x)
            Arblib.add!(ai, ai, tmp)

            # ai -= 2a0
            Arblib.mul_2exp!(tmp, Arblib.ref(res, 0), 1)
            Arblib.sub!(ai, ai, tmp)

            # ai /= 8
            Arblib.mul_2exp!(ai, ai, -3)

            res[2] = ai

            if deg >= 3
                # ai = (besselj(ν - 3, x) - besselj(ν + 3, x) - 6a1) / 48

                # ai = besselj(ν - 3, x)
                Arblib.sub!(tmp, ν, 3)
                Arblib.hypgeom_bessel_j!(ai, tmp, x)

                # tmp = besselj(ν + 3, x)
                # ai += tmp
                Arblib.add!(tmp, ν, 3)
                Arblib.hypgeom_bessel_j!(tmp, tmp, x)
                Arblib.sub!(ai, ai, tmp)

                # ai -= 6a1
                Arblib.mul!(tmp, Arblib.ref(res, 1), 6)
                Arblib.sub!(ai, ai, tmp)

                # ai /= 48
                Arblib.div!(ai, ai, 48)

                res[3] = ai

                if deg >= 4
                    x2 = x^2
                    x2mν2 = x2 - ν^2
                    for i = 4:deg
                        k = i - 4

                        # tmp = 2x * res[k + 1]
                        Arblib.mul_2exp!(tmp, x, 1)
                        Arblib.mul!(tmp, tmp, Arblib.ref(res, k + 1))

                        # ai = res[k] + tmp
                        Arblib.add!(ai, Arblib.ref(res, k), tmp)

                        # tmp = (k^2 + 4k + 4 - ν^2 + x^2)
                        # ai += tmp * res[k + 2]
                        Arblib.add!(tmp, x2mν2, k^2 + 4k + 4)
                        Arblib.addmul!(ai, tmp, Arblib.ref(res, k + 2))

                        # tmp = (2k^2 + 11k + 15) * x
                        # ai += tmp * res[k + 3]
                        Arblib.mul!(tmp, x, 2k^2 + 11k + 15)
                        Arblib.addmul!(ai, tmp, Arblib.ref(res, k + 3))

                        # tmp = x^2 * (k^2 + 7k + 12)
                        # ai /= tmp
                        Arblib.mul!(tmp, x2, k^2 + 7k + 12)
                        Arblib.div!(ai, ai, tmp)

                        # ai = -ai
                        Arblib.neg!(ai, ai)

                        res[i] = ai
                    end
                end
            end
        end
    end

    # Finally compose the Taylor series for the Bessel function with
    # that of z
    return compose_zero!(res, res, z)
end

SpecialFunctions.besselj0(z::ArbSeries) = besselj(zero(Arb), z)
SpecialFunctions.besselj1(z::ArbSeries) = besselj(one(Arb), z)
