using LambertWExp, LambertW, Test

"inverse of u = lambertwexp(x)"
lambexpinv(u) = log(u) + u


@testset "lambertwexp" begin


@test iszero(lambertwexp(-Inf))
@test lambertwexp(Inf) == Inf
@test isnan(lambertwexp(NaN))

@inferred lambertwexp(1)
@inferred lambertwexp(-1)


for x = -100 : 10 : 100
    @inferred lambertwexp(x)
    if iszero(x)
        @test lambexpinv(lambertwexp(x)) ≈ x  atol = 1e-15
    else
        @test lambexpinv(lambertwexp(x)) ≈ x
    end
end


for x = -10 : 1 : 10
    @test lambertwexp(x) ≈ lambertw(exp(x))
end


end  # testset
