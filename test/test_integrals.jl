
using DomainSets: ×

"Does `z` agree with `val` up to the given `threshold`?"
isaccurate(z, val, threshold = 1e-5) = abs(z-val) < threshold

function test_readme_examples()
    # Examples of README.md
    threshold = 1e-5
    @test isaccurate(integral(cos, 0..1.0), 0.8414709848078965)
    @test isaccurate(integral(exp(x) for x in 2..3), 12.69648082425702)
    @test isaccurate(integral(exp(x+y) for (x,y) in (0..1)^2), 2.9524924420120535)
    @test isaccurate(integrate(cos(x) for x in UnionDomain(0..1, 2..3))[1], 0.07329356604208204)
    @test isaccurate(integral( (sin(log(abs(t))) for t in  -1..1), LogSingPoint(0.0)), -1.0000000021051316)
    @test isaccurate(integral( ( exp(log(abs(x-y))) for (x,y) in (2..3) × (1..4) ), DiagonallySingular()),
2.333333333333333)
    @test isaccurate(integral(cos, ChebyshevTWeight()), 2.403939430634413)
    @test isaccurate(integral(cos(t)*1/sqrt(1-t^2) for t in  -1.0..1.0), 2.403939410869398)
    @test isaccurate(integrate(QuadAdaptive(atol=1e-3, rtol = 1e-3), t->cos(t^2), 0..10)[1], 0.6011251848111901, 1e-3)
    @test isaccurate(integral(Q_GaussLaguerre(10), cos), 0.5000005097999486)
    @test isaccurate(integral(cos(t)*exp(-t) for t in HalfLine()), 0.5)
end

function test_some_integrals()
    z1 = sin(1.0)
    @test abs(integral(cos, 0..1.0) - z1) < 1e-8
    I, E = integrate(cos, 0..1.0)
    @test abs(I-z1) < 1e-8
    @test E < 1e-8

    f2 = x -> exp(x[1]+x[2])
    z2 = hcubature(f2, (0,0), (1,1))[1]
    @test abs(integral(f2, (0..1)^2) - z2) < 1e-6
end

function test_integrals()
    test_some_integrals()
    test_readme_examples()
end
