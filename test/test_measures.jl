
function test_measures()
    m1 = LebesgueMeasure()
    @test domaintype(m1) == Float64
    @test DomainIntegrals.unsafe_weight(m1, 0.4) == 1

    m2 = UnitLebesgueMeasure()
    @test domaintype(m2) == Float64
    @test support(m2) == UnitInterval()
    @test weight(m2, 0.4) == 1
    @test weight(m2, 1.4) == 0

    m3 = LegendreMeasure()
    @test domaintype(m3) == Float64
    @test support(m3) == ChebyshevInterval()
    @test weight(m3, 0.4) == 1
    @test weight(m3, -0.4) == 1
    @test weight(m3, 1.4) == 0
    @test weight(m3, -1.4) == 0

    x = 0.5
    m4 = DiracMeasure(x)
    @test domaintype(m4) == Float64
    @test support(m4) == Point(x)
    @test point(m4) == x
    @test weight(m4, x) == Inf
    @test weight(m4, x+1) == 0
end
