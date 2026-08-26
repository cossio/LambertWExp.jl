# LambertWExp Julia Package

[![Build Status](https://travis-ci.org/cossio/LambertWExp.jl.svg?branch=master)](https://travis-ci.org/cossio/LambertWExp.jl)
[![Coverage Status](https://coveralls.io/repos/github/cossio/LambertWExp.jl/badge.svg?branch=master)](https://coveralls.io/github/cossio/LambertWExp.jl?branch=master)

Computes `W(e^x)` for real x and the principal branch of W, avoiding intermediate overflow.

Uses the fourth-order Fritsch iteration (Fritsch, Shafer & Crowell 1973; see also Veberic 2012, https://doi.org/10.1016/j.cpc.2012.07.008) adapted to the `W(e^x)` formulation, so that `e^x` is never formed explicitly.
