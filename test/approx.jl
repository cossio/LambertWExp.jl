using LambertW, ForwardDiff

import LambertWExp: lambertwexp_approx,
                    lambertwexp_approx_d1,
                    lambertwexp_approx_d2


@testset "lambertwexp_approx" begin  
    for x = -10:10
        exact = lambertwexp(x)

        @test lambertwexp_approx(x) ≈ exact  atol=0.3 rtol=0.2
        
        @test abs(lambertwexp_approx(x) - exact) ≤ 0.35
        @test abs(lambertwexp_approx(x) - exact) / exact ≤ 0.2
        
        @test lambertwexp_approx_d1(x) ≈ ForwardDiff.derivative(lambertwexp_approx, x)
        @test lambertwexp_approx_d2(x) ≈ ForwardDiff.derivative(lambertwexp_approx_d1, x)

        @test lambertwexp_approx_d1(x) > 0     # monotonic increasing
        @test lambertwexp_approx_d2(x) > 0     # convex
    end
end
