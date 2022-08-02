"""
    _compute_endpoints(a::Arf, b::Arf, x::Arb)

Given `a`, `b` and `x` with `x = Arb((a, b))` compute
```
c = a - midpoint(x)
d = b - midpoint(x)
```
rounded towards zero. It returns `c, c_exact, d, d_exact` where `c`
and `d` are as above and `c_exact` is true if `c` was computed exactly
and false if rounding was performed, similarly for `d_exact`.

This is used in [`extrema_series`](@ref), [`minimum_series`](@ref) and
[`maximum_series`](@ref) to compute the interval on which the extrema
of the polynomial should be enclosed.

If `c_exact` is false then an enclosure of the endpoint can be
computed with [`_enclose_inexact_endpoint`](@ref).
"""
function _compute_endpoints(a::Arf, b::Arf, x::Arb)
    c, d = zero(a), zero(b)
    xmid = Arblib.midref(x)
    c_exact = iszero(Arblib.sub!(c, a, xmid, rnd = Arblib.ArbRoundToZero))
    d_exact = iszero(Arblib.sub!(d, b, xmid, rnd = Arblib.ArbRoundToZero))
    return c, c_exact, d, d_exact
end

"""
    _enclose_inexact_endpoint(endpoint::Arf)

Compute a ball centered at `endpoint` with a radius of 1 ulp of
`endpoint`.

This is used to get an enclosure of the endpoint in case
[`_compute_endpoints`](@ref) gives an inexact result. The semantics of
Arb guarantees that the error of the endpoint is at most 1 ulp when
using directed rounding and a ball with `endpoint` as midpoint and a
radius set to the value of this 1 ulp will therefore always enclose
the exact result.
"""
function _enclose_inexact_endpoint(endpoint::Arf)
    res = Arb(endpoint)
    Arblib.set_ulp!(Arblib.radref(res), endpoint, prec = precision(endpoint))
    return res
end

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
    p = f(ArbSeries((x, 1), degree = degree + 1))

    # Check if f happens to be monotonic on the interval, in that case
    # evaluate on the endpoints only
    if !Arblib.contains_zero(Arblib.ref(p, 1))
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

    remainder = taylor_remainder(p, x)

    # If the remainder is finite the result will never be finite,
    # return the zeroth order enclosure
    if !isfinite(remainder)
        verbose && @info "non-finite remainder"
        res = maybe_abs(p[0])
        return res, copy(res), copy(res)
    end

    # Compute the Taylor series at the midpoint of x
    q = f(ArbSeries((Arblib.midref(x), 1); degree))

    # Compute interval to enclose extrema of Taylor series on
    c, c_exact, d, d_exact = _compute_endpoints(a, b, x)

    # Enclose the extrema of the Taylor series
    res = extrema_polynomial(q.poly, c, d; abs_value)

    if !c_exact
        # Compute value on an enclosure of the left endpoint
        y = maybe_abs(q(_enclose_inexact_endpoint(c)))
        Arblib.min!(res[1], res[1], y)
        Arblib.max!(res[2], res[2], y)
    end
    if !d_exact
        # Compute value on an enclosure of the right endpoint
        y = maybe_abs(q(_enclose_inexact_endpoint(d)))
        Arblib.min!(res[1], res[1], y)
        Arblib.max!(res[2], res[2], y)
    end

    return res[1] + remainder, res[2] + remainder, maybe_abs(q[0])
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
    p = f(ArbSeries((x, 1), degree = degree + 1))

    # Check if f happens to be monotonic on the interval, in that case
    # evaluate on the endpoints where the maximum could be attained
    if !Arblib.contains_zero(Arblib.ref(p, 1))
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
        elseif Arblib.ispositive(Arblib.ref(p, 1))
            # Minimum is attained at the left endpoint
            res = f(Arb(a))
        else
            # Minimum is attained at the right endpoint
            res = f(Arb(b))
        end
        return res, Arb(NaN, prec = precision(a))
    end

    remainder = taylor_remainder(p, x)

    # If the remainder is finite the result will never be finite,
    # return the zeroth order enclosure
    if !isfinite(remainder)
        verbose && @info "non-finite remainder"
        res = maybe_abs(p[0])
        return res, copy(res)
    end

    # Compute the Taylor series at the midpoint of x
    q = f(ArbSeries((Arblib.midref(x), 1); degree))

    # Compute interval to enclose minimum of Taylor series on
    c, c_exact, d, d_exact = _compute_endpoints(a, b, x)

    # Enclose the minimum of the Taylor series
    res = minimum_polynomial(q.poly, c, d; abs_value)

    if !c_exact
        # Compute value on an enclosure of the left endpoint
        Arblib.min!(res, res, maybe_abs(q(_enclose_inexact_endpoint(c))))
    end
    if !d_exact
        # Compute value on an enclosure of the right endpoint
        Arblib.min!(res, res, maybe_abs(q(_enclose_inexact_endpoint(d))))
    end

    return res + remainder, maybe_abs(q[0])
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
    p = f(ArbSeries((x, 1), degree = degree + 1))

    # Check if f happens to be monotonic on the interval, in that case
    # evaluate on the endpoints where the maximum could be attained
    if !Arblib.contains_zero(Arblib.ref(p, 1))
        verbose && @info "monotonic on interval - evaluate on endpoints"
        if abs_value
            # Maximum could be attained at either endpoint
            res = max(abs(f(Arb(a))), abs(f(Arb(b))))
        elseif Arblib.ispositive(Arblib.ref(p, 1))
            # Maximum is attained at the right endpoint
            res = f(Arb(b))
        else
            # Maximum is attained at the left endpoint
            res = f(Arb(a))
        end
        return res, Arb(NaN, prec = precision(a))
    end

    remainder = taylor_remainder(p, x)

    # If the remainder is finite the result will never be finite,
    # return the zeroth order enclosure
    if !isfinite(remainder)
        verbose && @info "non-finite remainder"
        res = maybe_abs(p[0])
        return res, copy(res)
    end

    # Compute the Taylor series at the midpoint of x
    q = f(ArbSeries((Arblib.midref(x), 1); degree))

    # Compute interval to enclose maximum of Taylor series on
    c, c_exact, d, d_exact = _compute_endpoints(a, b, x)

    # Enclose the maximum of the Taylor series
    res = maximum_polynomial(q.poly, c, d; abs_value)

    if !c_exact
        # Compute value on an enclosure of the left endpoint
        Arblib.max!(res, res, maybe_abs(q(_enclose_inexact_endpoint(c))))
    end
    if !d_exact
        # Compute value on an enclosure of the right endpoint
        Arblib.max!(res, res, maybe_abs(q(_enclose_inexact_endpoint(d))))
    end

    return res + remainder, maybe_abs(q[0])
end

"""
    enclosure_series(f, x::Arb; degree = 0, abs_value = false, verbose = false)

Convenience function for computing an enclosure of `f(x)` using
[`extrema_series`](@ref).

The enclosure is computed by finding the minimum and maximum value
using [`extrema_series`](@ref) and typically gives a tighter
enclosures than just evaluating `f(x)` directly.

It is equivalent to
```
Arb(ArbExtras.extrema_series(f, getinterval(x)...; degree, abs_value, verbose)[1:2])
```
but has the default value `degree = 0`, which is different than
[`extrema_series`](@ref). The reason for the default `degree = 0` is
that it can pick up on monotonicity of `f` and is therefore enough in
many cases.
"""
enclosure_series(f, x::Arb; degree::Integer = 0, abs_value = false, verbose = false) =
    Arb(extrema_series(f, getinterval(x)...; degree, abs_value, verbose)[1:2])
