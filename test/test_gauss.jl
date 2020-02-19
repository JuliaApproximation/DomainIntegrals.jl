
function test_gausslegendre()
    n = 4
    d = 2n
    qs = Q_GaussLegendre(n)
    @test integral(qs, x -> x^(d-1), 1..2) ≈ 2^d/d-1^d/d
    @test integral(BestRule(n), x -> x^(d-2), LegendreMeasure()) ≈ 2/(d-1)
end

function test_gaussjacobi()
    n = 4
    α = 1/2
    β = 3/4
    qs = Q_GaussJacobi(n, α, β)
    I = integral(x -> cos(x)*(1+x)^α*(1-x)^β, -1..1)
    @test abs(integral(qs, cos, -1..1) - I) < 1e-5
    @test abs(integral(BestRule(n), cos, JacobiMeasure(α, β)) -I) < 1e-5
end

function test_gausslaguerre()
    n = 15
    α = 1.2
    qs = Q_GaussLaguerre(n, α)
    I = integral(x->cos(x)*x^(α)*exp(-x), 0..10)
    @test abs(integral(qs, cos, HalfLine()) - I) < 1e-4
    @test abs(integral(BestRule(n), cos, LaguerreMeasure(α)) - I) < 1e-4
end

function test_gausshermite()
    n = 5
    qs = Q_GaussHermite(n)
    I = integral(x->cos(x)*exp(-x^2), -5..5)
    @test abs(integral(qs, cos, DomainSets.FullSpace{Float64}()) - I) < 1e-5
    @test abs(integral(BestRule(n), cos, HermiteMeasure()) - I) < 1e-5
end

function test_gauss()
    test_gausslegendre()
    test_gaussjacobi()
    test_gausslaguerre()
    test_gausshermite()
end
