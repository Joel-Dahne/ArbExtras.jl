"""
    _check_signs(pa, pb, values = ())

Check the signs of `pa`, `pb` and `values`. Returns `1` if all of them
definitely have the same sign, `-1` if there are at least two who
definitely have different signs and `0` otherwise.
"""
function _check_signs(pa, pb, values = ())
    has_negative = false
    has_positive = false
    has_zero = false

    for v in Base.Iterators.flatten(((pa, pb), values))
        sgn = Arblib.sgn_nonzero(v)

        has_negative |= sgn == -1
        has_positive |= sgn == 1
        has_zero |= sgn == 0

        has_negative && has_positive && return -1
    end

    if has_zero
        return 0
    else
        return 1
    end
end

"""
    extrema_polynomial(p::ArbPoly, a::Arf, b::Arf; abs_value = false, verbose = false)

Compute an enclosure of both the minimum and maximum of the polynomial
`p` on the interval `[a, b]` and return them as a 2-tuple.

The extrema are computed by finding the zeros of the derivative of `p`
using [`isolate_roots`](@ref) and refining them using
[`refine_root`](@ref), then evaluating `p` on the potential zeros as
well as the endpoints `a`, `b` of the interval.

If `abs_value = true` the compute the extrema of `abs(p(x))` on the
interval `[a, b]`. For the computation of the maximum this is mostly
the same, only difference is that we have to take the absolute value
of the evaluations. For the minimum we have to take into account that
the polynomial might cross zero, in which case the minimum is zero. To
see if `p` crosses zero we compare the sign of `p` at the endpoints
and all zeros of its derivative. If any of these have different signs
we are guaranteed that `p` crosses zero and the minimum is zero, if
all of them have the same sign we are instead guaranteed that the
polynomial doesn't cross zero.

If `verbose = true` then output information about the process.

TODO: Allow passing arguments to `isolate_roots`?
"""
function extrema_polynomial(p::ArbPoly, a::Arf, b::Arf; abs_value = false, verbose = false)
    isfinite(a) && isfinite(b) ||
        throw(ArgumentError("a and b must be finite, got a = $a and b = $b"))
    a <= b || throw(ArgumentError("must have a <= b, got a = $a and b = $b"))

    maybe_abs = abs_value ? abs : identity

    if a == b
        pa = maybe_abs(p(a))
        return pa, pa
    end

    # Compute values at endpoints
    pa, pb = p(a), p(b)

    # Flag to determine if the minimum value has already been
    # determined in the case when computing the absolute value.
    min_done = false

    # Handle the signs not being the same when computing the absolute
    # value
    if abs_value
        sgn = _check_signs(pa, pb)
        if sgn == -1
            # The sign differs - the minimum is zero
            verbose && @info "sign of endpoints differ - minimum is zero"

            min_value = zero(pa)
            min_done = true
        elseif sgn == 0
            # At least one of the enclosures contain zero, can't get a
            # better enclosure than the minimum of these two
            min_value = min(maybe_abs(pa), maybe_abs(pb))
            Arblib.nonnegative_part!(min_value, min_value)

            verbose &&
                @info "sign of endpoints undetermined - minimum on endpoints: $min_value"
            min_done = true
        end
    end

    # Compute the extrema on the endpoints
    min_endpoints, max_endpoints = minmax(maybe_abs(pa), maybe_abs(pb))

    verbose && @info "extrema on endpoints: $((min_endpoints, max_endpoints))"

    # Short circuit non-finite values and return NaN
    !isfinite(min_endpoints) ||
        !isfinite(max_endpoints) && return (
            Arb(NaN, prec = precision(min_endpoints)),
            Arb(NaN, prec = precision(min_endpoints)),
        )

    # Short circuit on constant or linear polynomial where extrema is
    # attained at endpoints (or possibly is zero if computing absolute
    # value)
    if Arblib.degree(p) <= 1
        if min_done
            return min_value, max_endpoints
        else
            return min_endpoints, max_endpoints
        end
    end

    # Compute the roots of the derivative
    dp = Arblib.derivative(p)
    roots, flags = isolate_roots(dp, a, b)
    roots = Arb.(roots) # We want the roots as balls rather than tuples

    verbose &&
        @info "found $(length(roots)) intervals with possible roots of which $(sum(flags)) unique"

    # Evaluate the polynomial on all the roots
    values = p.(roots)

    # Handle sign crossings when computing the absolute value
    if abs_value && !min_done && _check_signs(pa, pb, values) == -1
        verbose && @info "polynomial crosses zero - minimum is zero"
        min_value = zero(pa)
        min_done = true
    end

    if !min_done
        min_value = copy(min_endpoints)
    end
    max_value = copy(max_endpoints)
    for v in values
        v = maybe_abs(v)
        min_done || Arblib.min!(min_value, min_value, v)
        Arblib.max!(max_value, max_value, v)
    end

    abs_value && Arblib.nonnegative_part!(min_value, min_value)
    abs_value && Arblib.nonnegative_part!(max_value, max_value)

    verbose && @info "extrema with unrefined roots: $((min_value, max_value))"

    # Refine the unique roots where the extrema could occur
    count = 0
    for i in eachindex(roots)
        !(maybe_abs(values[i]) < max_value)

        if flags[i] && (
            !(min_done || maybe_abs(values[i]) > min_value) ||
            !(maybe_abs(values[i]) < max_value)
        )
            roots[i] = refine_root(dp, roots[i], strict = false)
            values[i] = p(roots[i])
            count += 1
        end
    end

    verbose && @info "refined $count roots"

    if abs_value && !min_done && _check_signs(pa, pb, values) == -1
        verbose && @info "polynomial crosses zero - minimum is zero"
        min_value = zero(pa)
        min_done = true
    end

    # Recompute extrema on the refined roots
    if !min_done
        min_value = copy(min_endpoints)
    end
    max_value = copy(max_endpoints)
    for v in values
        v = maybe_abs(v)
        min_done || Arblib.min!(min_value, min_value, v)
        Arblib.max!(max_value, max_value, v)
    end

    abs_value && Arblib.nonnegative_part!(min_value, min_value)
    abs_value && Arblib.nonnegative_part!(max_value, max_value)

    verbose && @info "extrema with refined roots: $((min_value, max_value))"

    return min_value, max_value
end

"""
    minimum_polynomial(p::ArbPoly, a::Arf, b::Arf; abs_value = false, verbose = false)

Compute an enclosure of the minimum of the polynomial `p` on the
interval `[a, b]`.

Takes the same arguments as [`extrema_polynomial`](@ref). The
algorithm is also the same except that it only looks for the minimum.
"""
function minimum_polynomial(p::ArbPoly, a::Arf, b::Arf; abs_value = false, verbose = false)
    isfinite(a) && isfinite(b) ||
        throw(ArgumentError("a and b must be finite, got a = $a and b = $b"))
    a <= b || throw(ArgumentError("must have a <= b, got a = $a and b = $b"))

    maybe_abs = abs_value ? abs : identity

    if a == b
        return maybe_abs(p(a))
    end

    # Compute values at endpoints
    pa, pb = p(a), p(b)

    # Handle the signs not being the same when computing the absolute
    # value
    if abs_value
        sgn = _check_signs(pa, pb)
        if sgn == -1
            # The sign differs return zero
            verbose && @info "sign of endpoints differ - minimum is zero"

            return zero(pa)
        elseif sgn == 0
            # At least one of the enclosures contain zero, return the
            # minimum of the two enclosures
            m = min(maybe_abs(pa), maybe_abs(pb))
            Arblib.nonnegative_part!(m, m)

            verbose && @info "sign of endpoints undetermined - minimum on endpoints: $m"

            return m
        end
    end

    # Compute the minimum on the endpoints
    m_endpoints = min(maybe_abs(pa), maybe_abs(pb))

    abs_value && Arblib.nonnegative_part!(m_endpoints, m_endpoints)

    verbose && @info "minimum on endpoints: $m_endpoints"

    # Short circuit non-finite values and return NaN
    !isfinite(m_endpoints) && return Arb(NaN, prec = precision(m_endpoints))

    # Short circuit on constant or linear polynomial where minimum is
    # attained at endpoints
    Arblib.degree(p) <= 1 && return m_endpoints

    # Compute the roots of the derivative
    dp = Arblib.derivative(p)
    roots, flags = isolate_roots(dp, a, b)
    roots = Arb.(roots) # We want the roots as balls rather than tuples

    verbose &&
        @info "found $(length(roots)) intervals with possible roots of which $(sum(flags)) unique"

    # Evaluate the polynomial on all the roots and find the minimum
    values = p.(roots)

    # Handle sign crossings when computing the absolute value
    if abs_value && _check_signs(pa, pb, values) == -1
        verbose && @info "polynomial crosses zero - minimum is zero"
        return zero(m_endpoints)
    end

    m = copy(m_endpoints)
    for v in values
        Arblib.min!(m, m, maybe_abs(v))
    end

    # If computing with absolute value keep only non-negative part
    abs_value && Arblib.nonnegative_part!(m, m)

    verbose && @info "minimum with unrefined roots: $m"

    # Refine the unique roots where the minimum could occur
    count = 0
    for i in eachindex(roots)
        if flags[i] && !(values[i] > m)
            roots[i] = refine_root(dp, roots[i], strict = false)
            values[i] = p(roots[i])
            count += 1
        end
    end

    verbose && @info "refined $count roots"

    # Handle sign crossings when computing the absolute value
    if abs_value && _check_signs(pa, pb, values) == -1
        verbose && @info "polynomial crosses zero - minimum is zero"
        return zero(m_endpoints)
    end

    # Recompute minimum on the refined roots
    m = copy(m_endpoints)
    for v in values
        Arblib.min!(m, m, maybe_abs(v))
    end

    # If computing with absolute value keep only non-negative part
    abs_value && Arblib.nonnegative_part!(m, m)

    verbose && @info "minimum with refined roots: $m"

    return m
end

"""
    maximum_polynomial(p::ArbPoly, a::Arf, b::Arf; abs_value = false, verbose = false)

Compute an enclosure of the maximum of the polynomial `p` on the
interval `[a, b]`.

Takes the same arguments as [`extrema_polynomial`](@ref). The
algorithm is also the same except that it only looks for the maximum.
"""
function maximum_polynomial(p::ArbPoly, a::Arf, b::Arf; abs_value = false, verbose = false)
    isfinite(a) && isfinite(b) ||
        throw(ArgumentError("a and b must be finite, got a = $a and b = $b"))
    a <= b || throw(ArgumentError("must have a <= b, got a = $a and b = $b"))

    maybe_abs = abs_value ? abs : identity

    if a == b
        return maybe_abs(p(a))
    end

    # Compute the maximum on the endpoints
    m_endpoints = max(maybe_abs(p(a)), maybe_abs(p(b)))

    # If computing with absolute value keep only non-negative part
    abs_value && Arblib.nonnegative_part!(m_endpoints, m_endpoints)

    verbose && @info "maximum on endpoints: $m_endpoints"

    # Short circuit non-finite values and return NaN
    !isfinite(m_endpoints) && return Arb(NaN, prec = precision(m_endpoints))

    # Short circuit on constant or linear polynomial where maximum is
    # attained at endpoints
    Arblib.degree(p) <= 1 && return m_endpoints

    # Compute the roots of the derivative
    dp = Arblib.derivative(p)
    roots, flags = isolate_roots(dp, a, b)
    roots = Arb.(roots) # We want the roots as balls rather than tuples

    verbose &&
        @info "found $(length(roots)) intervals with possible roots of which $(sum(flags)) unique"

    # Evaluate the polynomial on all the roots and find the maximum
    values = p.(roots)

    m = copy(m_endpoints)
    for v in values
        Arblib.max!(m, m, maybe_abs(v))
    end

    # If computing with absolute value keep only non-negative part
    abs_value && Arblib.nonnegative_part!(m, m)

    verbose && @info "maximum with unrefined roots: $m"

    # Refine the unique roots where the maximum could occur
    count = 0
    for i in eachindex(roots)
        if flags[i] && !(maybe_abs(values[i]) < m)
            roots[i] = refine_root(dp, roots[i], strict = false)
            values[i] = p(roots[i])
            count += 1
        end
    end

    verbose && @info "refined $count roots"

    # Recompute maximum on the refined roots
    m = copy(m_endpoints)
    for v in values
        Arblib.max!(m, m, maybe_abs(v))
    end

    # If computing with absolute value keep only non-negative part
    abs_value && Arblib.nonnegative_part!(m, m)

    verbose && @info "maximum with refined roots: $m"

    return m
end