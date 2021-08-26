"""
    bisect_interval(a::Arf, b::Arf; log_midpoint::Bool = false)

Returns two tuples, `(a, midpoint)` and `(midpoint, b)`, which
corresponds to splitting the interval in half.

If `log_midpoint = true` then the midpoint is computed in logarithmic
scale. If `a, b > 0` this sets the midpoint to `exp((log(a) + log(b))
/ 2)`. If `a, b < 0` the midpoint is `-exp((log(-a) + log(-b)) / 2)`.
If `a` or `b` is exactly zero then the midpoint is set to `exp(log(b)
/ 2)` and `-exp(log(-a) / 2)` respectively. If zero is contained in
the interior of the interval then currently a normal bisection is
performed (i.e. the midpoint is `(a + b) / 2`).

TODO: Think about how to handle a logarithmic bisection when the
interval contains zero

The value of `midpoint` is aliased in the two tuples and care should
therefore be taken if doing inplace operations on it.
"""
function bisect_interval(a::Arf, b::Arf; log_midpoint::Bool = false)
    if log_midpoint && a > 0
        c = Arb(a)
        d = Arb(b)
        Arblib.log!(c, c)
        Arblib.log!(d, d)
        Arblib.add!(c, c, d)
        Arblib.mul_2exp!(c, c, -1)
        Arblib.exp!(c, c)
        mid = midpoint(c)
    elseif log_midpoint && b < 0
        c = Arb(a)
        d = Arb(b)
        Arblib.neg!(c, c)
        Arblib.neg!(d, d)
        Arblib.log!(c, c)
        Arblib.log!(d, d)
        Arblib.add!(c, c, d)
        Arblib.mul_2exp!(c, c, -1)
        Arblib.exp!(c, c)
        Arblib.neg!(c, c)
        mid = midpoint(c)
    elseif log_midpoint && iszero(a) && iszero(b)
        mid = zero(a)
    elseif log_midpoint && iszero(a)
        d = Arb(b)
        Arblib.log!(d, d)
        Arblib.mul_2exp!(d, d, -1)
        Arblib.exp!(d, d)
        mid = midpoint(d)
    elseif log_midpoint && iszero(b)
        c = Arb(b)
        Arblib.neg!(c, c)
        Arblib.log!(c, c)
        Arblib.mul_2exp!(c, c, -1)
        Arblib.exp!(c, c)
        Arblib.neg!(c, c)
        mid = midpoint(c)
    else
        mid = a + b
        Arblib.mul_2exp!(mid, mid, -1)
    end

    return (a, mid), (mid, b)
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
"""
function check_tolerance(x::Arb; atol = nothing, rtol = nothing)
    isfinite(x) || return false
    isnothing(atol) && isnothing(rtol) && return true

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
