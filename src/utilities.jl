"""
    bisect_interval(a::Arf, b::Arf)

Returns two tuples, `(a, midpoint)` and `(midpoint, b)`, which
corresponds to splitting the interval in half.

The value of `midpoint` is aliased in the two tuples and care should
therefore be taken if doing inplace operations on it.
"""
function bisect_interval(a::Arf, b::Arf)
    midpoint = a + b
    Arblib.mul_2exp!(midpoint, midpoint, -1)
    return (a, midpoint), (midpoint, b)
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
check is skipped. If both of them are `nothing` then it always return
`true`
"""
function check_tolerance(x::Arb; atol = nothing, rtol = nothing)
    isnothing(atol) && isnothing(rtol) && return true
    !isfinite(x) && return false

    error = radius(Arb, x)
    Arblib.mul_2exp!(error, error, 1)

    !isnothing(atol) && error <= atol && return true

    if !isnothing(rtol)
        Arblib.contains_zero(x) && return iszero(error)

        return error <= rtol * abs(x)
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
end

"""
    check_interval(Bool, a::Arf, b::Arf)

Return `true` if `a` and `b` defines a finite interval `[a, b]`. It
checks that `a` and `b` are both finite and satisfy `a <= b`.
"""
check_interval(::Type{Bool}, a::Arf, b::Arf) = isfinite(a) && isfinite(b) && a <= b
