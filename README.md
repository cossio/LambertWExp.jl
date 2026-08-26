# LambertWExp Julia Package

[![Build Status](https://travis-ci.org/cossio/LambertWExp.jl.svg?branch=master)](https://travis-ci.org/cossio/LambertWExp.jl)
[![Coverage Status](https://coveralls.io/repos/github/cossio/LambertWExp.jl/badge.svg?branch=master)](https://coveralls.io/github/cossio/LambertWExp.jl?branch=master)

Computes `W(e^x)` for real x and the principal branch of W, avoiding intermediate overflow.

Based on https://github.com/jlapeyre/LambertW.jl.

## See also

`W(e^x)` is the [Wright omega function](https://en.wikipedia.org/wiki/Wright_omega_function) ω(x). The registered package [WrightOmega.jl](https://github.com/JuliaMath/WrightOmega.jl) implements it for real and complex arguments, with rules for several automatic differentiation packages. Consider using it instead of this package.

## TODO

Consider implementing Fritsch iteration, which should be faster (Veberic 2012).
