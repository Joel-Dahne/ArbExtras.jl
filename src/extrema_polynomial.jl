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
    _extrema_polynomial_low_degree(p::ArbPoly, a::Arf, b::Arf; abs_value = false)

Internal method used by [`extrema_polynomial`](@ref) for polynomials
of degree at most `2`.

It assumes that `check_interval(a, b)` holds but does not check it.

The method is split into separate cases depending on the degree of
`p`. It is optimized for performance and sacrifices readability.

**IMPROVE:** Add support for polynomials of degree `3`.
"""
function _extrema_polynomial_low_degree(p::ArbPoly, a::Arf, b::Arf; abs_value = false)
    deg = Arblib.degree(p)

    if deg <= 0
        # Constant polynomial
        res = p[0]

        if abs_value
            Arblib.abs!(res, res)
        end

        return res, copy(res)
    elseif deg == 1
        # Linear polynomial, extrema at endpoints
        p_a, p_b = p(a), p(b)

        if abs_value
            min_is_zero = Arblib.sgn_nonzero(p_a) * Arblib.sgn_nonzero(p_b) < 0

            Arblib.abs!(p_a, p_a)
            Arblib.abs!(p_b, p_b)

            if min_is_zero
                Arblib.max!(p_b, p_a, p_b)
                Arblib.zero!(p_a)
                return p_a, p_b
            else
                if p_a <= p_b
                    return p_a, p_b
                elseif p_a >= p_b
                    return p_b, p_a
                else
                    res_min = min(p_a, p_b)
                    Arblib.max!(p_b, p_a, p_b) # Reuse p_b
                    return res_min, p_b
                end
            end
        else
            sign = Arblib.sgn_nonzero(Arblib.ref(p, 1))

            if sign > 0
                # Increasing
                return p_a, p_b
            elseif sign < 0
                # Decreasing
                return p_b, p_a
            else
                res_min = min(p_a, p_b)
                Arblib.max!(p_b, p_a, p_b) # Reuse p_b
                return res_min, p_b
            end
        end
    elseif deg == 2
        # Compute value at endpoints
        p_a, p_b = p(a), p(b)

        # Compute value at unique critical point if it intersects
        # interval.
        p_root, root_intersects = let tmp = zero(p_a)
            if Arblib.contains_zero(Arblib.ref(p, 2))
                # Critical point could be anywhere or might not exist.
                # Just evaluate on the whole interval
                p_ab = Arb((a, b))
                Arblib.evaluate!(p_ab, p, p_ab)

                if abs_value
                    Arblib.abs!(p_a, p_a)
                    Arblib.abs!(p_b, p_b)
                    Arblib.abs!(p_ab, p_ab)
                end

                Arblib.min!(tmp, p_a, p_b)
                Arblib.min!(tmp, tmp, p_ab)
                Arblib.max!(p_b, p_a, p_b) # Reuse p_b
                Arblib.max!(p_b, p_b, p_ab)
                return tmp, p_b
            end

            # Compute root
            Arblib.div!(tmp, Arblib.ref(p, 1), Arblib.ref(p, 2))
            Arblib.neg!(tmp, tmp)
            Arblib.mul_2exp!(tmp, tmp, -1)

            # Check if root intersects the interval and in that case
            # evaluate, otherwise set to indeterminate.
            x = Arb((a, b))
            root_intersects = Arblib.overlaps(x, tmp)
            if root_intersects
                Arblib.intersection!(tmp, tmp, x)
                Arblib.evaluate!(tmp, p, tmp)
            end
            tmp, root_intersects
        end

        if abs_value
            # Check if polynomial has sign change
            p_a_sign = Arblib.sgn_nonzero(p_a)
            p_b_sign = Arblib.sgn_nonzero(p_b)
            p_root_sign = root_intersects ? Arblib.sgn_nonzero(p_root) : zero(p_a_sign)
            min_is_zero =
                any(<(0), (p_a_sign, p_b_sign, p_root_sign)) &&
                any(>(0), (p_a_sign, p_b_sign, p_root_sign))

            Arblib.abs!(p_a, p_a)
            Arblib.abs!(p_b, p_b)
            Arblib.abs!(p_root, p_root)

            if min_is_zero
                Arblib.max!(p_b, p_a, p_b) # Reuse p_b
                if root_intersects
                    Arblib.max!(p_b, p_b, p_root)
                end

                Arblib.zero!(p_a) # Reuse p_a

                return p_a, p_b
            else
                if root_intersects
                    res_min = min(p_a, p_b)
                    Arblib.min!(res_min, res_min, p_root)

                    Arblib.max!(p_b, p_a, p_b) # Reuse p_b
                    Arblib.max!(p_b, p_b, p_root)

                    return res_min, p_b
                else
                    Arblib.min!(p_root, p_a, p_b) # Reuse p_root
                    Arblib.max!(p_b, p_a, p_b)
                    return p_root, p_b
                end
            end
        else
            if root_intersects
                res_min = min(p_a, p_b)
                Arblib.min!(res_min, res_min, p_root)

                Arblib.max!(p_b, p_a, p_b) # Reuse p_b
                Arblib.max!(p_b, p_b, p_root)

                return res_min, p_b
            else
                Arblib.min!(p_root, p_a, p_b) # Reuse p_root
                Arblib.max!(p_b, p_a, p_b)
                return p_root, p_b
            end
        end
    else
        throw(ArgumentError("only supports polynomials of degree at most 2"))
    end
end


"""
    _minimum_polynomial_low_degree(p::ArbPoly, a::Arf, b::Arf; abs_value = false)

Internal method used by [`minimum_polynomial`](@ref) for polynomials
of degree at most `2`.

Similar to [`_extrema_polynomial_low_degree`](@ref) but only computes
minimum.
"""
function _minimum_polynomial_low_degree(p::ArbPoly, a::Arf, b::Arf; abs_value = false)
    deg = Arblib.degree(p)

    if deg <= 0
        # Constant polynomial
        res = p[0]

        if abs_value
            Arblib.abs!(res, res)
        end

        return res
    elseif deg == 1
        # Linear polynomial, extrema at endpoints
        if abs_value
            p_a, p_b = p(a), p(b)

            if Arblib.sgn_nonzero(p_a) * Arblib.sgn_nonzero(p_b) < 0
                return Arblib.zero!(p_a)
            end

            Arblib.abs!(p_a, p_a)
            Arblib.abs!(p_b, p_b)

            if p_a <= p_b
                return p_a
            elseif p_a >= p_b
                return p_b
            else
                return Arblib.min!(p_a, p_a, p_b) # Reuse p_a
            end
        else
            sign = Arblib.sgn_nonzero(Arblib.ref(p, 1))

            if sign > 0
                # Increasing
                return p(a)
            elseif sign < 0
                # Decreasing
                return p(b)
            else
                p_a, p_b = p(a), p(b)
                return Arblib.min!(p_a, p_a, p_b) # Reuse p_a
            end
        end
    elseif deg == 2
        # Compute value at endpoints
        p_a, p_b = p(a), p(b)

        if abs_value
            # If sign at endpoints differ minimum is zero
            p_a_sign = Arblib.sgn_nonzero(p_a)
            p_b_sign = Arblib.sgn_nonzero(p_b)

            if p_a_sign * p_b_sign < 0
                return Arblib.zero!(p_a)
            end
        end

        # Compute value at unique critical point if it intersects
        # interval.
        p_root, root_intersects = let tmp = zero(p_a)
            if Arblib.contains_zero(Arblib.ref(p, 2))
                # Critical point could be anywhere or might not exist.
                # Just evaluate on the whole interval
                p_ab = Arb((a, b))
                Arblib.evaluate!(p_ab, p, p_ab)

                if abs_value
                    Arblib.abs!(p_a, p_a)
                    Arblib.abs!(p_b, p_b)
                    Arblib.abs!(p_ab, p_ab)
                end

                Arblib.min!(tmp, p_a, p_b)
                Arblib.min!(tmp, tmp, p_ab)
                return tmp
            end

            # Compute root
            Arblib.div!(tmp, Arblib.ref(p, 1), Arblib.ref(p, 2))
            Arblib.neg!(tmp, tmp)
            Arblib.mul_2exp!(tmp, tmp, -1)

            # Check if root intersects the interval and in that case
            # evaluate, otherwise set to indeterminate.
            x = Arb((a, b))
            root_intersects = Arblib.overlaps(x, tmp)
            if root_intersects
                Arblib.intersection!(tmp, tmp, x)
                Arblib.evaluate!(tmp, p, tmp)
            end
            tmp, root_intersects
        end

        if abs_value
            # Check for sign change at critical point
            if root_intersects
                p_root_sign = Arblib.sgn_nonzero(p_root)
                if p_a_sign * p_root_sign < 0 || p_b_sign * p_root_sign < 0
                    return Arblib.zero!(p_a)
                end
            end

            Arblib.abs!(p_a, p_a)
            Arblib.abs!(p_b, p_b)
            Arblib.abs!(p_root, p_root)
        end

        Arblib.min!(p_a, p_a, p_b) # Reuse p_a

        if root_intersects
            Arblib.min!(p_a, p_a, p_root)
        end

        return p_a
    else
        throw(ArgumentError("only supports polynomials of degree at most 2"))
    end
end

"""
    _maximum_polynomial_low_degree(p::ArbPoly, a::Arf, b::Arf; abs_value = false)

Internal method used by [`minimum_polynomial`](@ref) for polynomials
of degree at most `2`.

Similar to [`_extrema_polynomial_low_degree`](@ref) but only computes
minimum.
"""
function _maximum_polynomial_low_degree(p::ArbPoly, a::Arf, b::Arf; abs_value = false)
    deg = Arblib.degree(p)

    if deg <= 0
        # Constant polynomial
        res = p[0]

        if abs_value
            Arblib.abs!(res, res)
        end

        return res
    elseif deg == 1
        # Linear polynomial, extrema at endpoints
        if abs_value
            p_a, p_b = p(a), p(b)

            Arblib.abs!(p_a, p_a)
            Arblib.abs!(p_b, p_b)

            if p_a <= p_b
                return p_b
            elseif p_a >= p_b
                return p_a
            else
                return Arblib.max!(p_a, p_a, p_b) # Reuse p_a
            end
        else
            sign = Arblib.sgn_nonzero(Arblib.ref(p, 1))

            if sign > 0
                # Increasing
                return p(b)
            elseif sign < 0
                # Decreasing
                return p(a)
            else
                p_a, p_b = p(a), p(b)
                return Arblib.max!(p_a, p_a, p_b) # Reuse p_a
            end
        end
    elseif deg == 2
        # Compute value at endpoints
        p_a, p_b = p(a), p(b)

        # Compute value at unique critical point if it intersects
        # interval.
        p_root, root_intersects = let tmp = zero(p_a)
            if Arblib.contains_zero(Arblib.ref(p, 2))
                # Critical point could be anywhere or might not exist.
                # Just evaluate on the whole interval
                p_ab = Arb((a, b))
                Arblib.evaluate!(p_ab, p, p_ab)

                if abs_value
                    Arblib.abs!(p_a, p_a)
                    Arblib.abs!(p_b, p_b)
                    Arblib.abs!(p_ab, p_ab)
                end

                Arblib.max!(tmp, p_a, p_b)
                Arblib.max!(tmp, tmp, p_ab)
                return tmp
            end

            # Compute root
            Arblib.div!(tmp, Arblib.ref(p, 1), Arblib.ref(p, 2))
            Arblib.neg!(tmp, tmp)
            Arblib.mul_2exp!(tmp, tmp, -1)

            # Check if root intersects the interval and in that case
            # evaluate, otherwise set to indeterminate.
            x = Arb((a, b))
            root_intersects = Arblib.overlaps(x, tmp)
            if root_intersects
                Arblib.intersection!(tmp, tmp, x)
                Arblib.evaluate!(tmp, p, tmp)
            end
            tmp, root_intersects
        end

        if abs_value
            Arblib.abs!(p_a, p_a)
            Arblib.abs!(p_b, p_b)
            Arblib.abs!(p_root, p_root)
        end

        Arblib.max!(p_a, p_a, p_b) # Reuse p_a

        if root_intersects
            Arblib.max!(p_a, p_a, p_root)
        end

        return p_a
    else
        throw(ArgumentError("only supports polynomials of degree at most 2"))
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

If `p` cannot be evaluated to a reasonable precision, because its
coefficients are too wide, it does in general not make sense to try
and isolate the roots. To check for this we compute `min(radius(p(a)),
radius(p(b)))` and `p(Arb((a, b)))`. If the quotient is larger than
`0.99`, meaning we don't even reduce the radius by one percent by
evaluating at the endpoints instead of the whole interval, we don't
try to isolate the roots but just evaluate the polynomial directly.

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
    check_interval(a, b)

    if a == b
        res = p(a)
        abs_value && Arblib.abs!(res, res)
        return res, copy(res)
    end

    Arblib.degree(p) <= 2 && return _extrema_polynomial_low_degree(p, a, b; abs_value)

    maybe_abs = abs_value ? abs : identity

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

    # Short circuit in case p cannot be evaluated to any meaningful
    # precision.
    pab = p(Arb((a, b)))
    if min(radius(pa), radius(pb)) / radius(pab) > Mag(0.99)
        verbose && @info "could not evaluate p to any meaningful precision"
        if min_done
            return min_value, max(max_endpoints, pab)
        else
            return min(min_endpoints, pab), max(max_endpoints, pab)
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
    check_interval(a, b)

    if a == b
        res = p(a)
        abs_value && Arblib.abs!(res, res)
        return res
    end

    Arblib.degree(p) <= 2 && return _minimum_polynomial_low_degree(p, a, b; abs_value)

    maybe_abs = abs_value ? abs : identity

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

    # Short circuit in case p cannot be evaluated to any meaningful
    # precision.
    pab = p(Arb((a, b)))
    if radius(m_endpoints) / radius(pab) > Mag(0.99)
        verbose && @info "could not evaluate p to any meaningful precision"
        return min(m_endpoints, pab)
    end

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
    check_interval(a, b)

    if a == b
        res = p(a)
        abs_value && Arblib.abs!(res, res)
        return res
    end

    Arblib.degree(p) <= 2 && return _maximum_polynomial_low_degree(p, a, b; abs_value)

    maybe_abs = abs_value ? abs : identity

    # Compute the maximum on the endpoints
    m_endpoints = max(maybe_abs(p(a)), maybe_abs(p(b)))

    # If computing with absolute value keep only non-negative part
    abs_value && Arblib.nonnegative_part!(m_endpoints, m_endpoints)

    verbose && @info "maximum on endpoints: $m_endpoints"

    # Short circuit non-finite values and return NaN
    !isfinite(m_endpoints) && return Arb(NaN, prec = precision(m_endpoints))

    # Short circuit in case p cannot be evaluated to any meaningful
    # precision.
    pab = p(Arb((a, b)))
    if radius(m_endpoints) / radius(pab) > Mag(0.99)
        verbose && @info "could not evaluate p to any meaningful precision"
        return max(m_endpoints, pab)
    end

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
