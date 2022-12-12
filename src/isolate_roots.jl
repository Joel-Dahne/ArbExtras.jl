"""
    _check_root_interval(f, a::Arf, b::Arf, fa_sign::Integer = Arblib.sgn_nonzero(f(Arb(a))), fb_sign::Integer = Arblib.sgn_nonzero(f(Arb(b))); check_unique = true, buffer1 = zero(Arb), buffer2 = zero(Arb))

Check if the function `f` has a zero on the interval `[a, b]`.

Returns `(maybe, unique)`. If `maybe` is `false` then `f` is proved to
not have a zero on the interval, if it's `true` then `f` might or
might not have a zero on the interval. If `unique` is `true` then `f`
is proved to have exactly one zero on the interval, if it's `false`
then no information is given.

The arguments `fa_sign` and `fb_sign` should be set to the sign of `f`
on the endpoints `a` and `b` according to how `Arblib.sgn_nonzero`
does it. They can be computed as
```
fa_sign = Arblib.sgn_nonzero(f(Arb(a)))
fb_sign = Arblib.sgn_nonzero(f(Arb(b)))
```

If `check_unique = false` then don't check for uniqueness, `unique`
will always be `false` in this case.

The function `f` should satisfy the same properties as for
[`isolate_roots`](@ref).

The arguments `buffer1` and `buffer2` can optionally be given to be
used as scratch space for the computations. This reduces the number of
allocations required.
"""
function _check_root_interval(
    f,
    a::Arf,
    b::Arf,
    fa_sign::Integer = Arblib.sgn_nonzero(f(Arb(a))),
    fb_sign::Integer = Arblib.sgn_nonzero(f(Arb(b)));
    check_unique::Bool = true,
    buffer1::Arb = zero(Arb),
    buffer2::Arb = zero(Arb), # Not used
)
    # If signs of endpoints could not be determined return immediately
    iszero(fa_sign) || iszero(fb_sign) && return true, false

    # buffer1 = Arb((a, b))
    Arblib.set_interval!(buffer1, a, b)

    maybe = Arblib.contains_zero(f(buffer1))

    maybe || return false, false

    check_unique || return true, false

    if fa_sign * fb_sign < 0
        df = Arblib.ref(f(ArbSeries((buffer1, 1))), 1)
        return true, !Arblib.contains_zero(df)
    else
        return true, false
    end
end

function _check_root_interval(
    p::ArbPoly,
    a::Arf,
    b::Arf,
    fa_sign::Integer = Arblib.sgn_nonzero(p(Arb(a))),
    fb_sign::Integer = Arblib.sgn_nonzero(p(Arb(b)));
    check_unique::Bool = true,
    buffer1::Arb = zero(Arb),
    buffer2::Arb = zero(Arb),
)
    # If signs of endpoints could not be determined return immediately
    iszero(fa_sign) || iszero(fb_sign) && return true, false

    # buffer1 = Arb((a, b))
    Arblib.set_interval!(buffer1, a, b)

    check_unique || return Arblib.contains_zero(p(buffer1)), false

    # buffer1, buffer2 = Arblib.evaluate2(p, Arb((a, b)))
    Arblib.evaluate2!(buffer1, buffer2, p, buffer1)

    Arblib.contains_zero(buffer1) || return false, false

    if fa_sign * fb_sign < 0
        unique = !Arblib.contains_zero(buffer2)
        return true, unique
    elseif fa_sign * fb_sign > 0
        maybe = Arblib.contains_zero(buffer2)
        return maybe, false
    else
        return true, false
    end
end

"""
    isolate_roots(f, a::Arf, b::Arf; depth = 10, check_unique = true, verbose = false)

Isolate the roots of `f` on the interval `[a, b]`.

Returns a tuple `(found, flags)` where `found` is a vector of
subintervals of `[a, b]` given by pairs of `Arf` and `flags` is a
vector of booleans. The output has the following properties:
- The function has no roots on interval outside of the subintervals in
  `found`.
- Subintervals are sorted in increasing order (with no overlap except
  possibly starting and ending with the same point).
- Subintervals with a flag `true` has exactly one (simple) root.
- Subintervals with a flag `false` may or may not contain roots.

If all flags are true then all roots of the function on interval have
been isolated. If there are output subintervals on which the existence
or nonexistence of roots could not be determined, the user may attempt
further searches on those subintervals (possibly with increased
precision and/or increased bounds for the breaking criteria). Note
that roots of multiplicity higher than one and roots located exactly
at endpoints cannot be isolated by the algorithm.

This method closely mimics the function `arb_calc_isolate_roots` in
Arb, both in behaviour and in implementation.

The function `f` should support evaluation on both `Arb` and
`ArbSeries` and should return an enclosure of the result in both
cases. The evaluation on `ArbSeries` is done to compute the derivative
to check uniqueness, this can be skipped by setting `check_unique =
false`, see below. Alternatively `f` can be an `ArbPoly`, in which
case some optimizations are done.

If `check_unique = false` then don't check if roots are unique. In
this case `flags` will always be false and it will simply continue to
split the subintervals until it reaches the maximum depth. The only
benefit to this is that the derivative of `f` doesn't have to be
computed so it can be used when the function of interest doesn't
support evaluating derivatives.

If `verbose = true` then output information after each iteration on
the number of found roots and the number of remaining intervals that
needs to be split further.

FUTURE WORK:
- Parallelize.
- Compute values on endpoints of subintervals once to not have to
  compute them multiple times.
- Should it use `BitVector` or `Vector{Bool}` for flags?
- The current implementation does a breadth-first-search. Consider
  rewriting it to do a depth-first-search instead to reduce memory
  usage in case of many splittings. This might however make it harder
  to parallelize.
"""
function isolate_roots(
    f,
    a::Arf,
    b::Arf;
    depth::Integer = 10,
    check_unique::Bool = true,
    verbose = false,
)
    check_interval(a, b)

    if a == b
        fa = f(Arb(a))
        if Arblib.contains_zero(fa)
            return [(a, b)], BitVector((iszero(fa),))
        else
            return Vector{NTuple{2,Arf}}(), BitVector()
        end
    end

    buffer1 = Arb(prec = Arblib._precision((a, b)))
    buffer2 = Arb(prec = Arblib._precision((a, b)))

    intervals = [(a, b)]

    Arblib.set!(buffer1, a)
    Arblib.set!(buffer2, b)
    if f isa ArbPoly
        fa_sign = Arblib.sgn_nonzero(Arblib.evaluate!(buffer1, f, buffer1))
        fb_sign = Arblib.sgn_nonzero(Arblib.evaluate!(buffer2, f, buffer2))
    else
        fa_sign = Arblib.sgn_nonzero(f(buffer1))
        fb_sign = Arblib.sgn_nonzero(f(buffer2))
    end
    sign_endpoints = [(fa_sign, fb_sign)]

    found = empty(intervals)
    flags = BitVector()

    iterations = 0
    while !isempty(intervals) && iterations < depth
        iterations += 1

        to_split = falses(length(intervals))
        for i in eachindex(intervals, sign_endpoints, to_split)
            maybe, unique = _check_root_interval(
                f,
                intervals[i]...,
                sign_endpoints[i]...;
                check_unique,
                buffer1,
                buffer2,
            )

            if unique
                push!(found, intervals[i])
                push!(flags, true)
            elseif maybe
                to_split[i] = true
            end
        end

        if iterations < depth
            intervals = bisect_intervals(intervals, to_split)

            # Compute signs of the endpoints for the newly bisected
            # intervals
            sign_endpoints_new = similar(intervals, eltype(sign_endpoints))
            i = 1
            for j in eachindex(to_split)
                if to_split[j]
                    mid = intervals[i][2]
                    Arblib.set!(buffer1, mid)
                    if f isa ArbPoly
                        Arblib.evaluate!(buffer1, f, buffer1)
                        fc_sign = Arblib.sgn_nonzero(buffer1)
                    else
                        fc_sign = Arblib.sgn_nonzero(f(buffer1))
                    end
                    sign_endpoints_new[i] = (sign_endpoints[j][1], fc_sign)
                    sign_endpoints_new[i+1] = (fc_sign, sign_endpoints[j][2])
                    i += 2
                end
            end
            sign_endpoints = sign_endpoints_new
        else
            # If we are on the last iteration don't split the intervals
            intervals = intervals[to_split]
        end

        verbose &&
            iterations < depth &&
            @info "isolate_roots iteration: $(lpad(iterations, 2)), " *
                  "found: $(lpad(length(found), 2)), " *
                  "remaining intervals: $(length(intervals) ÷ 2)"
    end

    verbose &&
        iterations == depth &&
        @info "isolate_roots iteration: $(lpad(iterations, 2)), " *
              "found: $(lpad(length(found), 2)), " *
              "remaining intervals: $(length(intervals))"

    if !isempty(intervals)
        append!(found, intervals)
        append!(flags, falses(length(intervals)))
    end

    p = sortperm(found, by = interval -> interval[1])

    return found[p], flags[p]
end
