"""
    _midpoint_interval(a::T, b::T) where {T<:Union{Arf,Arb}}

Return the midpoint of `a` and `b`.
"""
function _midpoint_interval(a::T, b::T) where {T<:Union{Arf,Arb}}
    mid = a + b
    Arblib.mul_2exp!(mid, mid, -1)

    return mid
end

"""
    _midpoint_interval_log!(buffer::Arb, a::Arf, b::Arf)

Return the logarithmic midpoint of `a` and `b` as described in
[`bisect_interval`](@ref).

The argument `buffer` is used as scratch space during the
computations.
"""
function _midpoint_interval_log!(buffer::Arb, a::Arf, b::Arf)
    cmpa = Arblib.cmp(a, 0)
    cmpb = Arblib.cmp(b, 0)

    if cmpa > 0 && cmpb > 0
        Arblib.set!(buffer, a)
        Arblib.mul!(buffer, buffer, b)

        Arblib.sqrt!(buffer, buffer)
        return midpoint(buffer)
    elseif cmpa < 0 && cmpb < 0
        Arblib.set!(buffer, a)
        Arblib.mul!(buffer, buffer, b)

        Arblib.sqrt!(buffer, buffer)
        Arblib.neg!(buffer, buffer)
        return midpoint(buffer)
    else
        mid = a + b
        Arblib.mul_2exp!(mid, mid, -1)

        return mid
    end
end

"""
    _midpoint_interval_log(a::Arf, b::Arf)

Return the logarithmic midpoint of `a` and `b` as described in
[`bisect_interval`](@ref).
"""
function _midpoint_interval_log(a::Arb, b::Arb)
    if Arblib.ispositive(a) && Arblib.ispositive(b)
        mid = a * b
        Arblib.sqrt!(mid, mid)

        return mid
    elseif Arblib.isnegative(a) && Arblib.isnegative(b)
        mid = a * b
        Arblib.sqrt!(mid, mid)
        Arblib.neg!(mid, mid)

        return mid
    else
        mid = a + b
        Arblib.mul_2exp!(mid, mid, -1)

        return mid
    end
end

"""
    bisect_interval(a::Arf, b::Arf; log_midpoint::Bool = false)

Compute the midpoint of `a and `b` and return two tuples, `(a,
midpoint)` and `(midpoint, b)`, which corresponds to splitting the
interval in half.

If `log_midpoint = true` then the midpoint is computed in logarithmic
scale. If `a, b > 0` this sets the midpoint to `exp((log(a) + log(b))
/ 2)`, can also be written as `sqrt(a * b)`. If `a, b < 0` we instead
get `-sqrt(a * b)`. If zero is contained in the interval then
currently a normal bisection is performed (i.e. the midpoint is `(a +
b) / 2`).

Logarithmic bisection can be useful if the function in consideration
has a logarithmic behaviour and the interval has numbers exponentially
close to zero, for example `[1e-10000, 1e-1]`. In that case normal
bisection could give very slow convergence.

TODO: Think about how to handle a logarithmic bisection when the
interval contains zero.

The value of `midpoint` is aliased in the two tuples and care should
therefore be taken if doing inplace operations on it.
"""
function bisect_interval(a::T, b::T; log_midpoint::Bool = false) where {T<:Union{Arf,Arb}}
    if log_midpoint
        if T == Arf
            mid = _midpoint_interval_log!(Arb(prec = Arblib._precision((a, b))), a, b)
        else
            mid = _midpoint_interval_log(a, b)
        end
    else
        mid = _midpoint_interval(a, b)
    end

    return (a, mid), (mid, b)
end

"""
    bisect_interval_recursive(a::T, b::T, depth::Integer; log_midpoint::Bool = false) where {T<:Union{Arf,Arb}}

Recursively bisect the interval `[a, b]` a number of `depth` times.
The bisection is done using [`bisect_interval`](@ref) and it returns a
vector with `2^depth` intervals.

It accepts the same arguments as [`bisect_interval`](@ref).
"""
function bisect_interval_recursive(
    a::T,
    b::T,
    depth::Integer;
    log_midpoint::Bool = false,
) where {T<:Union{Arf,Arb}}
    depth >= 0 || Throw(ArgumentError("depth needs to be non-negative, got $depth"))

    res = Vector{NTuple{2,T}}(undef, 2^depth)
    res[1] = (a, b)
    @inbounds for i = 1:depth
        for j in reverse(1:2^(i-1))
            res[2j-1], res[2j] = bisect_interval(res[j]...; log_midpoint)
        end
    end

    return res
end

"""
    bisect_intervals(intervals::Vector{NTuple{2,T}}, to_bisect::Union{BitVector,Vector{Bool}}; log_midpoint = false) where {T<:Union{Arf,Arb}}

Bisect all intervals in `intervals` for which `to_bisect` is true and
return a vector with the bisected intervals.

See also [`bisect_interval`](@ref).
"""
function bisect_intervals(
    intervals::Vector{NTuple{2,T}},
    to_bisect::Union{BitVector,Vector{Bool}};
    log_midpoint::Bool = false,
) where {T<:Union{Arf,Arb}}
    length(intervals) == length(to_bisect) || throw(
        ArgumentError("intervals and to_bisect should have the same number of elements"),
    )

    res = Vector{eltype(intervals)}(undef, 2sum(to_bisect))

    isempty(res) && return res

    if T == Arf && log_midpoint
        buffer = Arb()
    end

    index = 1
    for i in eachindex(to_bisect)
        if to_bisect[i]
            a, b = intervals[i]

            if log_midpoint
                if T == Arf
                    mid = _midpoint_interval_log!(buffer, a, b)
                else
                    mid = _midpoint_interval_log(a, b)
                end
            else
                mid = _midpoint_interval(a, b)
            end

            res[index], res[index+1] = (a, mid), (mid, b)

            index += 2
        end
    end

    return res
end

"""
    check_tolerance(x::Arb; atol = nothing, rtol = nothing)

Return `true` if `x` satisfies the absolute tolerance `atol` or the
relative tolerance `rtol`.

The absolute tolerance is satisfied if the diameter of `x` is less
than or equal to `atol`.

The relative tolerance is satisfied if the diameter of `x` divided by
the absolute value of `x` is less than or equal to `rtol`. If `x`
contains zero this is only satisfied if `x` is exactly zero.

Both `atol` and `rtol` can be given as `nothing`, in which case that
check is skipped. If both of them are `nothing` then it returns `true`
if `x` is finite.

It always returns true if `x` is finite and the radius is zero.
"""
function check_tolerance(x::Arb; atol = nothing, rtol = nothing)
    isfinite(x) || return false
    isnothing(atol) && isnothing(rtol) && return true

    iszero(Arblib.radref(x)) && return true

    # Radius is always non-zero from here

    error = radius(Arb, x)
    Arblib.mul_2exp!(error, error, 1)

    !isnothing(atol) && !iszero(atol) && error <= atol && return true

    if !isnothing(rtol) && !iszero(rtol) && !Arblib.contains_zero(x)
        bound = rtol * x
        return error <= Arblib.abs!(bound, bound)
    else
        return false
    end
end

"""
    check_interval(a::Arf, b::Arf)

Throw an error if `a` and `b` does not define a finite interval `[a,
b]`. It checks that `a` and `b` are both finite and satisfy `a <= b`.
"""
function check_interval(a::Arf, b::Arf)
    isfinite(a) && isfinite(b) ||
        throw(ArgumentError("a and b must be finite, got a = $a and b = $b"))
    a <= b || throw(ArgumentError("must have a <= b, got a = $a and b = $b"))
    return nothing
end

"""
    check_interval(Bool, a::Arf, b::Arf)

Return `true` if `a` and `b` defines a finite interval `[a, b]`. It
checks that `a` and `b` are both finite and satisfy `a <= b`.
"""
check_interval(::Type{Bool}, a::Arf, b::Arf) = isfinite(a) && isfinite(b) && a <= b

"""
    format_interval(a, b; digits = 5)

Return a nicely formatted string representing the interval `[a, b]`.

It prints the interval either on the form `[a, b]` or on the form `[c
+/- d]` depending on the width of the width of the interval.

If `a` or `b` is non-finite it always prints an interval of the form
`[a, b]`. Otherwise it checks if the relative error is less than
`10^-digits`, if it is then it prints it as a ball, otherwise as an
interval. When printed in interval form it prints at most `digits`
digits.
"""
function format_interval(a::Arf, b::Arf; digits = 5)
    if isfinite(a) && isfinite(b)
        ball = Arb((a, b))
        if check_tolerance(ball, rtol = Arb(10)^(-digits))
            return string(ball)
        end
    end

    return "[" *
           Base.MPFR._string(convert(BigFloat, a), digits) *
           ", " *
           Base.MPFR._string(convert(BigFloat, b), digits) *
           "]"
end

"""
    taylor_remainder(p::ArbSeries, x::Arb)

Compute `p[degree] * (x - midpoint(x))^degree`, where `degree` is the
degree of `p`

This corresponds to the Lagrange form of the remainder term in a
Taylor expansion around the point `midpoint(x)` of degree `degree -
1` valid on the full interval `x`.

For example, for a function `f`, a ball `x`, any `degree >= 0` and
```
p = f(ArbSeries((x, 1), degree = degree + 1))
q = f(ArbSeries((midpoint(x), 1), degree = degree))
```
we have for any `y ∈ x` that
```
f(y) ∈ q(y) + taylor_remainder(p, x)
```
"""
function taylor_remainder(p::ArbSeries, x::Arb)
    restterm = zero(x)
    Arblib.set!(Arblib.radref(restterm), Arblib.radref(x))
    Arblib.pow!(restterm, restterm, unsigned(Arblib.degree(p)))
    Arblib.mul!(restterm, restterm, Arblib.ref(p, Arblib.degree(p)))
    return restterm
end
