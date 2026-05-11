"""
    _taylor_with_lagrange_remainder(f, x0::T, interval::T; degree::Integer, enclosure_degree = 0) where {T<:Union{Arb,Acb}}

Compute the Taylor expansion of `f` at the point `x0` of degree
`degree` with the last term being the remainder term in Lagrange form,
which ensures that the truncated version gives an enclosure of `f(x)`
for all `x ∈ interval`.

We require that `x0 ∈ interval`.

For wide values of `interval`, it computes a tighter enclosure of the
remainder term using [`ArbExtras.enclosure_series`](@ref). The degree
used for this can be set with `enclosure_degree`. Setting it to a
negative number computes the remainder directly instead. For `Acb`,
`enclosure_degree` is accepted but has no effect; the remainder term
is always computed directly.
"""
function _taylor_with_lagrange_remainder(
    f,
    x0::T,
    interval::T;
    degree::Integer,
    enclosure_degree::Integer = 0,
) where {T<:Union{Arb,Acb}}
    Arblib.contains(interval, x0) || throw(
        ArgumentError(
            "expected x0 to be contained in interval, got x0 = $x0, interval = $interval",
        ),
    )

    TSeries = T == Arb ? ArbSeries : AcbSeries

    if x0 == interval
        # In this case, we can compute the full expansion directly
        return f(TSeries((x0, 1); degree))
    end

    # Compute the expansion without the remainder term
    res = f(TSeries((x0, 1), degree = degree - 1))

    # Make room for the remainder term
    res = TSeries(res; degree)

    # Compute the remainder term
    if T == Arb && enclosure_degree >= 0
        # We compute a tighter enclosure with the help of enclosure_series
        res[degree] =
            enclosure_series(
                derivative_function(f, degree),
                interval,
                degree = enclosure_degree,
            ) / factorial(degree)
    else
        res[degree] = f(TSeries((interval, 1); degree))[degree]
    end

    return res
end

"""
    fx_div_x(f, x::Union{Arb,Acb}, order::Integer = 1; extra_degree = 0, enclosure_degree = 0, force = false)
    fx_div_x(f, x::Union{ArbSeries,AcbSeries}, order::Integer = 1; extra_degree = 0, enclosure_degree = 0, force = false)

Compute an enclosure of `f(x) / x^order` for a function `f` with a
zero of the given order at the origin. This can be used to compute
enclosures for functions with removable singularities.

Setting `extra_degree` to a value higher than `0` makes it use a
higher order expansion to enclose the value. This can give tighter
bounds in many cases. The argument `enclosure_degree` is passed to
[`_taylor_with_lagrange_remainder`](@ref).

If `force = false`, it requires that the enclosure of `f` at zero is
exactly zero. If `f` is known to be exactly zero at zero but the
enclosure might be wider, it can be forced to be zero by setting `force
= true`.

This function is based on Lemma A.1 in Appendix A of

> Dahne, J., & Gómez-Serrano, J. (2023). Highest Cusped Waves for the
> Burgers–Hilbert Equation. *Archive for Rational Mechanics and
> Analysis*, 247(5). https://doi.org/10.1007/s00205-023-01904-6
"""
function fx_div_x(
    f,
    x::Union{Arb,Acb},
    order::Integer = 1;
    extra_degree::Integer = 0,
    enclosure_degree::Integer = 0,
    force::Bool = false,
)
    Arblib.contains_zero(x) || throw(ArgumentError("x must contain zero, got x = $x"))
    order >= 1 || throw(ArgumentError("order must be at least 1, got $order"))
    extra_degree >= 0 ||
        throw(ArgumentError("extra_degree must be non-negative, got $extra_degree"))

    if iszero(x)
        # If x is exactly zero, there is no need to compute extra terms
        # in the expansion
        extra_degree = 0
    end

    expansion = _taylor_with_lagrange_remainder(
        f,
        zero(x),
        x,
        degree = order + extra_degree;
        enclosure_degree,
    )

    isfinite(expansion) || return Arblib.indeterminate!(zero(x))

    if force
        for i = 0:(order-1)
            Arblib.contains_zero(Arblib.ref(expansion, i)) || error(
                "coefficient $i of the expansion does not contain zero, got $(expansion[i])",
            )
            expansion[i] = 0
        end
    end

    return Arblib.shift_right!(
        typeof(expansion)(
            degree = Arblib.degree(expansion) - order,
            prec = precision(expansion),
        ),
        expansion,
        order,
    )(
        x,
    )
end

function fx_div_x(
    f,
    x::T,
    order::Integer = 1;
    extra_degree::Integer = 0,
    enclosure_degree::Integer = 0,
    force::Bool = false,
) where {T<:Union{ArbSeries,AcbSeries}}
    x0 = x[0]
    Arblib.contains_zero(x0) ||
        throw(ArgumentError("x[0] must contain zero, got x[0] = $x0"))
    order >= 1 || throw(ArgumentError("order must be at least 1, got $order"))
    extra_degree >= 0 ||
        throw(ArgumentError("extra_degree must be non-negative, got $extra_degree"))

    expansion = _taylor_with_lagrange_remainder(
        f,
        zero(x0),
        x0,
        degree = Arblib.degree(x) + order + extra_degree;
        enclosure_degree,
    )

    if !isfinite(expansion)
        indeterminate_res = zero(x)
        ind_Arb = Arblib.indeterminate!(zero(x0))
        for i = 0:Arblib.degree(x)
            indeterminate_res[i] = ind_Arb
        end
        return indeterminate_res
    end

    if force
        for i = 0:(order-1)
            Arblib.contains_zero(Arblib.ref(expansion, i)) || error(
                "coefficient $i of the expansion does not contain zero, got $(expansion[i])",
            )
            expansion[i] = 0
        end
    end

    expansion_div_x = Arblib.shift_right!(
        T(degree = Arblib.degree(expansion) - order, prec = precision(expansion)),
        expansion,
        order,
    )

    # Set the result to the Taylor series of f(x) / x^order at x0
    res = T(degree = Arblib.degree(x))
    for i = 0:Arblib.degree(res)
        res[i] = expansion_div_x(x0) / factorial(i)
        Arblib.derivative!(expansion_div_x, expansion_div_x)
    end

    return compose_zero(res, x)
end
