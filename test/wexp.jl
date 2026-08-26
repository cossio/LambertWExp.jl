using LambertW, ForwardDiff


"inverse of u = lambertwexp(x)"
lambexpinv(u) = log(u) + u


@testset "lambertwexp, special values and type stability" begin
    @test iszero(lambertwexp(-Inf))
    @test lambertwexp(Inf) == Inf
    @test isnan(lambertwexp(NaN))
    @test isnan(lambertwexp(Float16(NaN)))

    for x in -100 : 100
        @inferred lambertwexp(x)
        @inferred lambertwexp(Float16(x))
        @inferred lambertwexp(Float64(x))
    end

    for x in [NaN, Inf, -Inf]
        @inferred lambertwexp(x)
        @inferred lambertwexp(Float16(x))
    end
end


@testset "lambertwexp numerical value" begin
    for x = -100 : 1 : 10000
        @inferred lambertwexp(x)
        @test lambexpinv(lambertwexp(x)) ≈ x  atol = 1e-15
    end

    for x in [1e4, 1e5, 1e6]
        @test lambexpinv(lambertwexp(x)) ≈ x
    end
end


@testset "lambertwexp accuracy against BigFloat reference" begin
    #= W(e^x) computed exactly via LambertW at 256-bit precision,
    covering the small-argument series region (x ≤ -2), the series
    about zero (-2 < x ≤ 1), and the asymptotic region (x > 1),
    including x > 709 where e^x overflows Float64. =#
    #= x between about -745 and -708 is avoided: there e^x is subnormal,
    so a 1-ulp difference in exp exceeds any reasonable rtol. =#
    xs = [-1000; -708; -300; -50:0.25:8; 10 .^ (1:15); 700; 710]
    setprecision(BigFloat, 256) do
        for x in Float64.(xs)
            ref = Float64(lambertw(exp(big(x))))
            @test lambertwexp(x) ≈ ref  rtol = 4eps(Float64)
        end
    end
end


@testset "lambertwexp huge arguments" begin
    #= For w ≳ 9.5e153 the Fritsch q-factor overflows Float64; these
    values exercise the overflow-robust form of the update. e^x is far
    beyond BigFloat range here, so accuracy is checked through the
    defining equation w + log(w) = x instead of a lambertw reference. =#
    setprecision(BigFloat, 256) do
        for x in [9.4e153, 9.5e153, 1e154, 1e200, 1e308, prevfloat(Inf)]
            w = lambertwexp(x)
            @test isfinite(w)
            @test big(w) + log(big(w)) ≈ big(x)  rtol = 4eps(Float64)
        end
    end
end


@testset "lambertwexp derivatives" begin
    for w = 1 : 1000
        x = w + log(w)
        @test lambertwexp_d1_from_W(w) ≈ ForwardDiff.derivative(lambertwexp, x)
        @test lambertwexp_d2_from_W(w) ≈ ForwardDiff.derivative(lambertwexp_d1_from_W ∘ lambertwexp, x)
    end
end