# ArbExtras.jl

This package extends [Arblib](https://github.com/kalmarek/Arblib.jl)
with some extra functionality. It's mainly intended for personal use
in my research and should not be seen as a stable Julia package. It
will likely never be added to the Julia repository. It can be seen as
a follow up of [ArbTools](https://github.com/Joel-Dahne/ArbTools.jl)
but based on Arblib instead of
[Nemo](https://github.com/wbhart/Nemo.jl) for the Arb interface.

## Installation
The package is not in the general Julia repository but can be
installed through the package manager with
``` julia
pkg> add https://github.com/Joel-Dahne/ArbExtras.jl
```

To see if it (some) things work correctly you can run the tests with
``` julia
pkg> test ArbExtras
```
