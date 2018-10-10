using LambertW, ForwardDiff


"inverse of u = lambertwexp(x)"
lambexpinv(u) = log(u) + u


@testset "lambertwexp, special values and type stability" begin
    @test iszero(lambertwexp(-Inf))
    @test lambertwexp(Inf) == Inf
    @test isnan(lambertwexp(NaN))

    @inferred lambertwexp(1)
    @inferred lambertwexp(-1)

    for x in [1e4, 1e5, 1e6]
        @test lambexpinv(lambertwexp(x)) ≈ x
    end
end


@testset "lambertwexp numerical value" begin
    for x = -100 : 10 : 100
        @inferred lambertwexp(x)
        if iszero(x)
            @test lambexpinv(lambertwexp(x)) ≈ x  atol = 1e-15
        else
            @test lambexpinv(lambertwexp(x)) ≈ x
        end
    end

    for x = -10 : 1 : 10
        @test lambexpinv(lambertwexp(x)) ≈ x  atol = 1e-15
        @test lambertwexp(x) ≈ lambertw(exp(x))
    end
end


@testset "lambertwexp derivatives" begin
    for w = 1 : 10
        x = w + log(w)
        @test lambertwexp_d1_from_W(w) ≈ ForwardDiff.derivative(lambertwexp, x)
        @test lambertwexp_d2_from_W(w) ≈ ForwardDiff.derivative(lambertwexp_d1_from_W ∘ lambertwexp, x)
    end
end