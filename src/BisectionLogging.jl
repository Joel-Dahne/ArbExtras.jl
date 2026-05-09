"""
    BisectionLogging{T}

Data structure for logging information in an algorithm applying an
adaptive bisection approach.

It consists of two vectors. One vector with subintervals represented
by pairs of `Arf` values and one vector with data of type `T`
associated to each subinterval. If `logger::BisectionLogging` is the
logger, then the subintervals are given by `logger.intervals` and the
associated data by `logger.data`.

The subintervals are disjoint, sorted and cover the entire interval on
which the bisection was applied.

The data associated to each interval depends on the method from which
this is returned. For e.g. [`minimum_enclosure`](@ref) and
[`maximum_enclosure`](@ref) the data is of type `Arb` and represents
an enclosure of the minimum, respectively maximum, on each subinterval.

# Extended help

## Implementing a method using `BisectionLogging`

It is the responsibility of the method logging the data to ensure that
the above specifications are satisfied. During the logging process it
is usually convenient to keep the subintervals unsorted and only sort
them at the end. Data can then be added to the `logger` with e.g.

```
push!(logger.intervals, intervals[i])
push!(logger.data, values[i])
```

Before the logger is returned, it should be sorted so that the
subintervals are ordered. This can be done with

```
sort_logger!(logger)
```
"""
struct BisectionLogging{T}
    intervals::Vector{NTuple{2,Arf}}
    data::Vector{T}
end

BisectionLogging{T}() where {T} = BisectionLogging(NTuple{2,Arf}[], T[])

function sort_logger!(logger::BisectionLogging)
    perm = sortperm(logger.intervals, by = first)
    permute!(logger.intervals, perm)
    permute!(logger.data, perm)
    return logger
end
