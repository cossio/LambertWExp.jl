# LambertWExp Julia Package

[![CI](https://github.com/cossio/LambertWExp.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/cossio/LambertWExp.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/cossio/LambertWExp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cossio/LambertWExp.jl)

Computes `W(e^x)` for real x and the principal branch of W, avoiding intermediate overflow.

Uses the fourth-order Fritsch iteration (Fritsch, Shafer & Crowell 1973; see also Veberic 2012, https://doi.org/10.1016/j.cpc.2012.07.008) adapted to the `W(e^x)` formulation, so that `e^x` is never formed explicitly.

## Comparison with WrightOmega.jl

`W(e^x)` is the [Wright omega function](https://en.wikipedia.org/wiki/Wright_omega_function) ω(x). The registered package [WrightOmega.jl](https://github.com/JuliaMath/WrightOmega.jl) computes the same function, using Fukushima's (2020) piecewise minimax rational approximations for real arguments (no iteration at all) and TOMS Algorithm 917 for complex arguments, and ships differentiation rules for ChainRulesCore, ForwardDiff, Enzyme, Mooncake and Symbolics.

The numbers below compare both packages against a 256-bit `BigFloat` reference solving `w + log(w) = x` (Julia 1.12.1, x86-64, WrightOmega v0.1.0; 2000–4000 log-spaced samples per region).

**Accuracy**, as error in ulps (max / mean):

| region | LambertWExp | WrightOmega |
|---|---|---|
| x ∈ [−745, −708] (ω subnormal) | 0.5 / 0.3 | flushes to `0.0` |
| x ∈ [−708, −2] | 1.6 / 0.35 | 3.7 / 1.2 |
| x ∈ [−2, 2] | 3.4 / 0.3 | 3.9 / 0.75 |
| x ∈ [2, 700] | 0.6 / 0.25 | 3.9 / 0.9 |
| x ∈ [700, 10³⁰⁸] (`e^x` overflows) | 0.5 / 0.01 | 2.5 / 0.03 |

**Speed**: WrightOmega costs about one `exp` call (~5–10 ns) for any real argument, since it only evaluates one rational function. The Fritsch iteration here costs 2–20× that (~10–130 ns), worst near x = 0 where the series seed is least accurate and more iterations are needed.

**Behavioral differences**:

- WrightOmega flushes results below `floatmin(Float64)` to exactly `0.0`, so it returns `0.0` for all x < −708.4. This package returns correctly rounded subnormals down to `5.0e-324` at x = −745.
- This package iterates in the argument's own precision, so `BigFloat` arguments get full-precision results; WrightOmega routes every real type through its `Float64` kernel, so `BigFloat` arguments only get `Float64` accuracy.
- WrightOmega supports complex arguments and integrates with automatic differentiation via package extensions; this package is real-only and instead exports the manual derivative helpers `lambertwexp_d1_from_W` and `lambertwexp_d2_from_W`.

In short: prefer WrightOmega.jl for speed, complex arguments, or AD integration; this package is useful when you need subnormal results in the deep tail or arbitrary-precision arguments.
