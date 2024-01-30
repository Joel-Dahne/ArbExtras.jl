# This file contains temporary methods which will likely be added to
# Arblib later on

"""
    iscpx(p::Union{ArbPoly,AcbPoly,ArbSeries,AcbSeries})

Return `true` if `p` is of the form `c + x`.
"""
iscpx(p::Union{ArbPoly,AcbPoly,ArbSeries,AcbSeries}) =
    length(Arblib.cstruct(p)) == 2 && isone(Arblib.ref(p, 1))

"""
    compose_zero!(res::T, p::T, q::T) where {T<:Union{ArbPoly,AcbPoly,ArbSeries,AcbSeries}}

In-place version of [`compose_zero`](@ref).
"""
function compose_zero!(res::T, p::T, q::T) where {T<:Union{ArbPoly,AcbPoly}}
    iscpx(q) && return Arblib.set!(res, p)

    if !iszero(q) && !iszero(Arblib.ref(q, 0))
        q = copy(q)
        q[0] = 0
    end

    return Arblib.compose!(res, p, q)
end

function compose_zero!(res::T, p::T, q::T) where {T<:Union{ArbSeries,AcbSeries}}
    iscpx(q) && return Arblib.set!(res, p)

    if !iszero(q) && !iszero(Arblib.ref(q, 0))
        q = copy(q)
        q[0] = 0
    end

    return Arblib.compose_series!(res, p, q, length(res))
end

"""
    compose_zero(p::T, q::T) where {T<:Union{ArbPoly,AcbPoly,ArbSeries,AcbSeries}}

Compute `compose(p, q)`, but with the constant term of `q` set to zero.
"""
function compose_zero(p::T, q::T) where {T<:Union{ArbPoly,AcbPoly}}
    res = T(prec = Arblib._precision(p, q))

    return compose_zero!(res, p, q)
end

function compose_zero(p::T, q::T) where {T<:Union{ArbSeries,AcbSeries}}
    res = T(degree = Arblib._degree(p, q), prec = Arblib._precision(p, q))

    return compose_zero!(res, p, q)
end
