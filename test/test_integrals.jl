
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
    @test isaccurate(integral( ( exp(log(abs(x-y))) for (x,y) in (2..3) × (1..4) ), SingularDiagonal()),
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
    @test abs(integral(cos,uniondomain(0..1.0, 2..3.0)) - (sin(1)+sin(3)-sin(2))) < 1e-6

    f2 = x -> exp(x[1]+x[2])
    z2 = hcubature(f2, (0,0), (1,1))[1]
    @test abs(integral(f2, (0..1)^2) - z2) < 1e-6

    @test integral(exp, Point(0.5)) == 0.0

    # a mapped domain
    @test abs(integral(cos, DomainSets.MappedDomain(LinearMap(1/2), 0..1))-sin(2)) < 1e-6

    # diagonal singularity
    @test abs(integral((log(abs(x-y)) for (x,y) in UnitSquare()), SingularDiagonal())+1.5) < 1e-5

    # a mapped weight
    testmap_d = 2.0..4.0
    testmap_μ = lebesguemeasure(testmap_d)
    testmap_m = LinearMap(4.0)
    I1 = integral(cos, testmap_d, testmap_μ)
    @test abs(integral(t -> cos(inverse(testmap_m,t)), testmap_m.(2.0..4.0),
        DomainIntegrals.mappedmeasure(testmap_m, testmap_μ)) - I1)  < 1e-6
    @test abs(integral(t -> cos(inverse(testmap_m,t)), MappedDomain(inverse(testmap_m), 2.0..4.0),
        DomainIntegrals.mappedmeasure(testmap_m, testmap_μ)) - I1)  < 1e-6

    @test integrate(cos, EmptySpace()) === (0.0, 0.0)
    @test integrate(DomainIntegrals.Integrand{ComplexF64}(cos), EmptySpace()) === (0.0+0.0im, 0.0)
    @inferred integrate(cos, EmptySpace())

    # product weights
    μ1 = ChebyshevTWeight()
    @test abs(integral(x->cos(x[1]+x[2]), (-0.5..0.5)^2, productmeasure(μ1, μ1)) - 1.0049593511549983) < 1e-6
    @test abs(integral(x->cos(x[1]+x[2]), (-0.5..0.5)^3, productmeasure(μ1, μ1, μ1)) - 1.0523909717390105) < 1e-6

    # parametric domains
    @test abs(integral(x->cos(x[1]+x[2]), UnitCircle()) - 3.51314344095576) < 1e-6
    @test abs(integral(x->cos(x[1]+x[2]), 5*UnitCircle()) - 9.41394510) < 1e-6

    g(z) = 5/z
    @test abs(1/(2*pi*im)*integral(g, ComplexUnitCircle()) - 5) < 1e-6

    strategy = QuadAdaptive(maxevals=1e6)
    @test abs(integral(strategy, x->1, UnitDisk()) - pi) < 1e-5
end

function test_integrals()
    test_some_integrals()
    test_readme_examples()
end
