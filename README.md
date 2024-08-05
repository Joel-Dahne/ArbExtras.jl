# ArbExtras.jl
This package extends [Arblib](https://github.com/kalmarek/Arblib.jl)
with some methods for enclosing roots, enclosing extrema and computing
integrals.

The package started development during my PhD with the goal to give
reusable implementations of algorithms that were used in my research.
The type of methods that are implemented and which type of input they
are optimized for has been heavily influenced by the projects that I
worked on during that time.

## Functionality
The main functionality provided by the package are methods for
isolating roots and enclosing extrema of univariate real functions. In
addition to this is provides a rudimentary integral routine and a few
minor things.

### Enclosing roots
The package provides three different functions for isolating and
enclosing roots of univariate, real valued functions:

- `isolate_roots(f, a::Arf, b::Arf)` is used for isolating all roots
  of a function `f` on the interval `[a, b]`.
- `refine_root(f, root::Arb)` is used for refining an initial
  enclosure `root` of a zero to `f`. It makes use of the interval
  Newton method.
- `refine_root_bisection(f, a, b)` is used for refining an initial
  enclosure `[a, b]` of a root of `f`. It requires that `f(a)` and
  `f(b)` have different signs. It makes use if bisection, checking the
  sign of the midpoint in each iteration.

See the documentation of the respective methods for more details, in
particular the available keyword arguments. The first two methods
require access to the derivative of `f`, this is computed using
`ArbSeries` and for that reason `f` needs to support evaluation of
both `Arb` and `ArbSeries`. The `refine_root_bisection` doesn't
require access to the derivative and can therefore be used when this
is not available.

### Enclosing extrema
The package provides three families of methods:

- `extrema_polynomial`, `minimum_polynomial` and `maximum_polynomial`
- `extrema_series`, `minimum_series` and `maximum_series`
- `extrema_enclosure`, `minimum_enclosure` and `maximum_enclosure`

As indicated by the names they compute either the minimum, the maximum
or both.

- `extrema_polynomial(p::ArbPoly, a::Arf, b::Arb)` is used for
  enclosing the extrema of a polynomial `p` on the interval `[a, b]`.
  This is done by enclosing the roots of the derivative of `p`.
- `extrema_series(f, a::Arf, b::Arb; degree::Integer = 8)` is used for
  enclosing the extrema of a function `f` on the interval `[a, b]`.
  This is done by computing a Taylor series of the given degree. The
  extrema of the series is computed using `extrema_polynomial` and
  then an enclosure of the remainder term is added.
- `extrema_enclosure(f, a::Arf, b::Arb)` is also
  used for enclosing the extrema of a function `f` on the interval
  `[a, b]`. In this case it is done by iteratively bisecting the
  interval and using `extrema_series` on each subinterval. The
  bisection is continued until the result satisfies the required
  tolerance.

See the documentation of the respective methods for more details, in
particular the available keyword arguments. The `minimum_` and
`maximum_` versions have the same interface, the only difference being
that they return either the minimum or maximum instead of both.

For `enclosure_series` it is also possible to call it like
`enclosure_series(f, x::Arb)`, in which case it computes an enclosure
of `f(x)` by computing the extrema and taking the union of them.

The package also provides `bounded_by(f, a::Arf, b::Arf, C::Arf)`
which checks if the function `f` is bounded by `C` on the interval
`[a, b]`. It is similar to `maximum_enclosure` but only bisects enough
to be able to check that it is bounded by `C`. See also the
`ubound_tol` and `lbound_tol` arguments to `extrema_enclosure`.

### Other functionality
In addition to the above the package also provides some smaller
things:

- `integrate(f, a::Arb, b::Arb)` can be used to compute the integral
  of `f` from `a` to `b` using a Gauss-Legendre quadrature of order 2,
  combined with bisection. In general the integrator from Arb,
  `Arblib.integrate` performs much better and should be preferred.
  This method does however have the benefit that it does not require
  evaluation on complex balls, instead it needs access to the fourth
  derivative of `f`, which is computed using `ArbSeries`. So it can be
  of use if there is no implementation of `f(::Acb)`, but there is one
  for `f(::ArbSeries)`.
- `SpecialFunctions.besselj(ν::Arb, z::ArbSeries)` and
  `SpecialFunctions.bessely(ν::Arb, z::ArbSeries)` for computing
  Taylor expansions of the Bessel functions.
- Utility functions which are mainly for internal use, but could
  be useful in other cases as well , see their documentation for
  details.
  - `bisect_interval`, `bisect_interval_recursive` and `bisect_intervals`
  - `check_tolerance`
  - `check_interval`
  - `format_interval`
  - `taylor_remainder`
  - `enclosure_ubound`, `enclosure_lbound` and `enclosure_getinterval`
  - `derivative_function`
- Some "temporary" functions that should possibly be added to
  Arblib.jl or even the Arb library itself in the future
  - `iscpx`
  - `compose_zero` and `compose_zero!`
