
function test_integrals()
    z1 = sin(1.0)
    @test abs(integral(cos, 0..1.0) - z1) < 1e-8
    I, E = quadrature(cos, 0..1.0)
    @test abs(I-z1) < 1e-8
    @test E < 1e-8

    f2 = x -> exp(x[1]+x[2])
    z2 = hcubature(f2, (0,0), (1,1))[1]
    @test abs(integral(f2, (0..1)^2) - z2) < 1e-6
end
