
function test_integrals()
    z1 = sin(1.0)
    @test abs(integral(cos, 0..1.0) - z1) < 1e-8
    I, E = quadrature(cos, 0..1.0)
    @test abs(I-z1) < 1e-8
    @test E < 1e-8
end
