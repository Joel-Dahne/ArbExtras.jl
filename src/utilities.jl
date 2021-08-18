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
