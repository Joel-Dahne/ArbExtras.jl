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

    # We want to enclose the extrema between c and d. However due to
    # rounding the values of c and d might not be known exactly. We
    # pass the smaller interval to extrema_polynomial, i.e. the one
    # given by taking the upper bound of c and the lower bound of d.
    # If c or d are not exact we also have to include the missing
    # parts of the endpoints.
    # TODO: In some extremely rare cases we might have ubound(c) >
    # lbound(d), in this case an exception will be thrown.
    c = a - midpoint(Arb, x)
    d = b - midpoint(Arb, x)

    # Enclose the extrema of the Taylor series
    res = extrema_polynomial(ArbPoly(q), ubound(c), lbound(d); abs_value)

    if !iszero(Arblib.radref(c))
        # Enclose extrema on [c, ubound(c)]
        y = maybe_abs(q(union(c, ubound(Arb, c))))
        Arblib.min!(res[1], res[1], y)
        Arblib.max!(res[2], res[2], y)
    end
    if !iszero(Arblib.radref(d))
        # Enclose extrema on [lbound(d), d]
        y = maybe_abs(q(union(d, lbound(Arb, d))))
        Arblib.min!(res[1], res[1], y)
        Arblib.max!(res[2], res[2], y)
    end

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
    # evaluate on the endpoints where the maximum could be attained
    if !Arblib.contains_zero(p[1])
        verbose && @info "monotonic on interval - evaluate on endpoints"
        if abs_value
            fa, fb = f(Arb(a)), f(Arb(b))
            sgn = _check_signs(fa, fb)
            if sgn == -1
                # The sign differs return zero
                verbose && @info "sign of endpoints differ - minimum is zero"

                res = zero(fa)
            else
                res = min(abs(fa), abs(fb))
                Arblib.nonnegative_part!(res, res)
            end
        elseif Arblib.ispositive(p[1])
            # Minimum is attained at the left endpoint
            res = f(Arb(a))
        else
            # Minimum is attained at the right endpoint
            res = f(Arb(b))
        end
        return res, Arb(NaN, prec = precision(a))
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

    # We want to enclose the minimum between c and d. However due to
    # rounding the values of c and d might not be known exactly. We
    # pass the smaller interval to minimum_polynomial, i.e. the one
    # given by taking the upper bound of c and the lower bound of d.
    # If c or d are not exact we also have to include the missing
    # parts of the endpoints.
    # FIXME: In some extremely rare cases we might have ubound(c) >
    # lbound(d), in this case an exception will be thrown.
    c = a - midpoint(Arb, x)
    d = b - midpoint(Arb, x)

    # Enclose the minimum of the Taylor series
    res = minimum_polynomial(ArbPoly(q), ubound(c), lbound(d); abs_value)

    if !iszero(Arblib.radref(c))
        # Enclose maximum on [c, ubound(c)]
        Arblib.min!(res, res, maybe_abs(q(union(c, ubound(Arb, c)))))
    end
    if !iszero(Arblib.radref(d))
        # Enclose maximum on [lbound(d), d]
        Arblib.min!(res, res, maybe_abs(q(union(d, lbound(Arb, d)))))
    end

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
    # evaluate on the endpoints where the maximum could be attained
    if !Arblib.contains_zero(p[1])
        verbose && @info "monotonic on interval - evaluate on endpoints"
        if abs_value
            # Maximum could be attained at either endpoint
            res = max(abs(f(Arb(a))), abs(f(Arb(b))))
        elseif Arblib.ispositive(p[1])
            # Maximum is attained at the right endpoint
            res = f(Arb(b))
        else
            # Maximum is attained at the left endpoint
            res = f(Arb(a))
        end
        return res, Arb(NaN, prec = precision(a))
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

    # We want to enclose the maximum between c and d. However due to
    # rounding the values of c and d might not be known exactly. We
    # pass the smaller interval to maximum_polynomial, i.e. the one
    # given by taking the upper bound of c and the lower bound of d.
    # If c or d are not exact we also have to include the missing
    # parts of the endpoints.
    # FIXME: In some extremely rare cases we might have ubound(c) >
    # lbound(d), in this case an exception will be thrown.
    c = a - midpoint(Arb, x)
    d = b - midpoint(Arb, x)

    # Enclose the maximum of the Taylor series
    res = maximum_polynomial(ArbPoly(q), ubound(c), lbound(d); abs_value)

    if !iszero(Arblib.radref(c))
        # Enclose maximum on [c, ubound(c)]
        Arblib.max!(res, res, maybe_abs(q(union(c, ubound(Arb, c)))))
    end
    if !iszero(Arblib.radref(d))
        # Enclose maximum on [lbound(d), d]
        Arblib.max!(res, res, maybe_abs(q(union(d, lbound(Arb, d)))))
    end

    return res + restterm, maybe_abs(q[0])
end
