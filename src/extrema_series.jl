"""
    extrema_series(f, a::Arf, b::Arf; degree, abs_value, verbose)

Compute both the minimum and maximum of the function `f` on the
interval `[a, b]`. It returns a three tuple `(fmin, fmax, fmid)` where
`fmin` and `fmax` are enclosures of the minimum and maximum
respectively and `fmid` is an enclosure of `f` evaluated on the
midpoint of the interval.

The extrema are computed by computing the extrema of the Taylor series
of `f` on the midpoint of the interval and then adding an error term.
The degree of the Taylor series used is set with the `degree`
argument.

The reason that `fmid` is returned is that as part of computing the
Taylor series we get this value for free. In some cases the Taylor
series for at the midpoint is however not computed, in this case
`fmid` will just be `NaN` instead.

If `abs_value = true` then compute the extrema of `abs(f(x))` on the
interval `[a, b]`. For the computation of the maximum this is mostly
the same, only difference is that we have to compute the maximum of
the absolute value of the Taylor series. For the minimum we have to
take into account that the polynomial might cross zero, in which case
the minimum is zero.

If `verbose = true` then output information about the process.

This method is mainly intended for internal use in the function
[`extrema_enclosure`](@ref) but can be used directly as well.

TODO: Handle crossing zero when computing absolute minimum better.
Currently we always add the rest term in the end.
"""
function extrema_series(
    f,
    a::Arf,
    b::Arf;
    degree::Integer = 8,
    abs_value = false,
    verbose = false,
)
    check_interval(a, b)

    maybe_abs = abs_value ? abs : identity

    if a == b
        fa = maybe_abs(f(Arb(a)))
        return fa, fa, fa
    end

    x = Arb((a, b))

    # Compute rest term of the Taylor series
    p = f(ArbSeries([x, one(x)], degree = degree + 1))

    # Check if f happens to be monotonic on the interval, in that case
    # evaluate on the endpoints only
    if !Arblib.contains_zero(p[1])
        verbose && @info "monotonic on interval - evaluate on endpoints"
        fa, fb = f(Arb(a)), f(Arb(b))
        if abs_value
            sgn = _check_signs(fa, fb)
            if sgn == -1
                # The sign differs return zero
                verbose && @info "sign of endpoints differ - minimum is zero"

                return zero(fa), max(abs(fa), abs(fb)), Arb(NaN, prec = precision(a))
            else
                res = minmax(abs(fa), abs(fb))
                Arblib.nonnegative_part!(res[1], res[1])

                return res..., Arb(NaN, prec = precision(a))
            end
        else
            return min(fa, fb), max(fa, fb), Arb(NaN, prec = precision(a))
        end
    end

    restterm = let y = zero(x)
        Arblib.set!(Arblib.radref(y), Arblib.radref(x))
        y^(degree + 1) * p[degree+1]
    end

    # If the restterm is finite the result will never be finite,
    # return the zeroth order enclosure
    if !isfinite(restterm)
        verbose && @info "non-finite restterm"
        res = maybe_abs(p[0])
        return res, copy(res), copy(res)
    end

    # Compute the Taylor series at the midpoint of x
    q = f(ArbSeries([Arblib.midref(x), one(x)]; degree))

    # FIXME: The computation of the endpoints here needs to be exact!
    # This can always be done by computing at sufficiently high
    # precision but is not currently done
    if a - Arblib.midref(x) != Arb(a) - Arblib.midpoint(Arb, x) ||
       b - Arblib.midref(x) != Arb(b) - Arblib.midpoint(Arb, x)
        throw(AssertionError("endpoints needs to be evaluated at higher precision"))
    end

    # Enclose the maximum of the Taylor series
    res = extrema_polynomial(
        ArbPoly(q),
        a - Arblib.midref(x),
        b - Arblib.midref(x);
        abs_value,
    )

    return res[1] + restterm, res[2] + restterm, maybe_abs(q[0])
end
"""
    minimum_series(f, a::Arf, b::Arf; degree, abs_value, verbose)

Compute the minimum of the function `f` on the interval `[a, b]`.

Takes the same arguments as [`extrema_series`](@ref). The algorithm is
also the same except that it only looks for the minimum.

This method is mainly intended for internal use in the function
[`minimum_enclosure`](@ref) but can be used directly as well.
"""
function minimum_series(
    f,
    a::Arf,
    b::Arf;
    degree::Integer = 8,
    abs_value = false,
    verbose = false,
)
    check_interval(a, b)

    maybe_abs = abs_value ? abs : identity
    if a == b
        fa = maybe_abs(f(Arb(a)))
        return fa, fa
    end

    x = Arb((a, b))

    # Compute rest term of the Taylor series
    p = f(ArbSeries([x, one(x)], degree = degree + 1))

    # Check if f happens to be monotonic on the interval, in that case
    # evaluate on the endpoints only
    if !Arblib.contains_zero(p[1])
        verbose && @info "monotonic on interval - evaluate on endpoints"
        fa, fb = f(Arb(a)), f(Arb(b))

        if abs_value
            sgn = _check_signs(fa, fb)
            if sgn == -1
                # The sign differs return zero
                verbose && @info "sign of endpoints differ - minimum is zero"

                return zero(fa), Arb(NaN, prec = precision(a))
            else
                m = min(abs(fa), abs(fb))
                Arblib.nonnegative_part!(m, m)

                return m, Arb(NaN, prec = precision(a))
            end
        else
            return min(fa, fb), Arb(NaN, prec = precision(a))
        end
    end

    restterm = let y = zero(x)
        Arblib.set!(Arblib.radref(y), Arblib.radref(x))
        y^(degree + 1) * p[degree+1]
    end

    # If the restterm is finite the result will never be finite,
    # return the zeroth order enclosure
    if !isfinite(restterm)
        verbose && @info "non-finite restterm"
        res = maybe_abs(p[0])
        return res, copy(res)
    end

    # Compute the Taylor series at the midpoint of x
    q = f(ArbSeries([Arblib.midref(x), one(x)]; degree))

    # FIXME: The computation of the endpoints here needs to be exact!
    # This can always be done by computing at sufficiently high
    # precision but is not currently done
    if !iszero(Arblib.radref(Arb(a) - Arblib.midpoint(Arb, x))) ||
       !iszero(Arblib.radref(Arb(b) - Arblib.midpoint(Arb, x)))
        throw(AssertionError("endpoints needs to be evaluated at higher precision"))
    end

    # Enclose the maximum of the Taylor series
    res = minimum_polynomial(
        ArbPoly(q),
        a - Arblib.midref(x),
        b - Arblib.midref(x);
        abs_value,
    )

    return res + restterm, maybe_abs(q[0])
end

"""
    maximum_series(f, a::Arf, b::Arf; degree, abs_value, verbose)

Compute the maximum of the function `f` on the interval `[a, b]`.

Takes the same arguments as [`extrema_series`](@ref). The algorithm is
also the same except that it only looks for the maximum.

This method is mainly intended for internal use in the function
[`maximum_enclosure`](@ref) but can be used directly as well.
"""
function maximum_series(
    f,
    a::Arf,
    b::Arf;
    degree::Integer = 8,
    abs_value = false,
    verbose = false,
)
    check_interval(a, b)

    maybe_abs = abs_value ? abs : identity

    if a == b
        fa = maybe_abs(f(Arb(a)))
        return fa, fa
    end

    x = Arb((a, b))

    # Compute rest term of the Taylor series
    p = f(ArbSeries([x, one(x)], degree = degree + 1))

    # Check if f happens to be monotonic on the interval, in that case
    # evaluate on the endpoints only
    if !Arblib.contains_zero(p[1])
        verbose && @info "monotonic on interval - evaluate on endpoints"
        fa, fb = f(Arb(a)), f(Arb(b))
        return max(maybe_abs(fa), maybe_abs(fb)), Arb(NaN, prec = precision(a))
    end

    restterm = let y = zero(x)
        Arblib.set!(Arblib.radref(y), Arblib.radref(x))
        y^(degree + 1) * p[degree+1]
    end

    # If the restterm is finite the result will never be finite,
    # return the zeroth order enclosure
    if !isfinite(restterm)
        verbose && @info "non-finite restterm"
        res = maybe_abs(p[0])
        return res, copy(res)
    end

    # Compute the Taylor series at the midpoint of x
    q = f(ArbSeries([Arblib.midref(x), one(x)]; degree))

    # FIXME: The computation of the endpoints here needs to be exact!
    # This can always be done by computing at sufficiently high
    # precision but is not currently done
    if a - Arblib.midref(x) != Arb(a) - Arblib.midpoint(Arb, x) ||
       b - Arblib.midref(x) != Arb(b) - Arblib.midpoint(Arb, x)
        throw(AssertionError("endpoints needs to be evaluated at higher precision"))
    end

    # Enclose the maximum of the Taylor series
    res = maximum_polynomial(
        ArbPoly(q),
        a - Arblib.midref(x),
        b - Arblib.midref(x);
        abs_value,
    )

    return res + restterm, maybe_abs(q[0])
end
