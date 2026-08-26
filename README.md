# LambertWExp Julia Package

[![CI](https://github.com/cossio/LambertWExp.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/cossio/LambertWExp.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/cossio/LambertWExp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cossio/LambertWExp.jl)

Computes `W(e^x)` for real x and the principal branch of W, avoiding intermediate overflow.

Based on https://github.com/jlapeyre/LambertW.jl.

## See also

`W(e^x)` is the [Wright omega function](https://en.wikipedia.org/wiki/Wright_omega_function) ω(x). The registered package [WrightOmega.jl](https://github.com/JuliaMath/WrightOmega.jl) implements it for real and complex arguments, with rules for several automatic differentiation packages. Consider using it instead of this package.

## TODO

Consider implementing Fritsch iteration, which should be faster (Veberic 2012).
