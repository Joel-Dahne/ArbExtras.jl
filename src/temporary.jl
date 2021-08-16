# This file contains temporary methods which will likely be added to
# Arblib later on

# Without this the method is ambiguous
Base.:^(x::ArbSeries, y::Rational) = x^Arb(y, prec = precision(x))
