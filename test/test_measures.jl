
function test_measures()
    m1 = LebesgueMeasure()
    @test !isdiscrete(m1)
    @test iscontinuous(m1)
    @test !isnormalized(m1)
    @test domaintype(m1) == Float64
    @test DomainIntegrals.unsafe_weight(m1, 0.4) == 1

    m2 = UnitLebesgueMeasure()
    @test !isdiscrete(m2)
    @test iscontinuous(m2)
    @test !isnormalized(m2)
    @test domaintype(m2) == Float64
    @test support(m2) == UnitInterval()
    @test weight(m2, 0.4) == 1
    @test weight(m2, 1.4) == 0

    m3 = LegendreMeasure()
    @test !isdiscrete(m3)
    @test iscontinuous(m3)
    @test !isnormalized(m3)
    @test domaintype(m3) == Float64
    @test support(m3) == ChebyshevInterval()
    @test weight(m3, 0.4) == 1
    @test weight(m3, -0.4) == 1
    @test weight(m3, 1.4) == 0
    @test weight(m3, -1.4) == 0

    x = 0.5
    m4 = DiracMeasure(x)
    @test isdiscrete(m4)
    @test !iscontinuous(m4)
    @test isnormalized(m4)
    @test domaintype(m4) == Float64
    @test support(m4) == Point(x)
    @test point(m4) == x
    @test weight(m4, x) == Inf
    @test weight(m4, x+1) == 0

    m5 = GaussianMeasure{SVector{2,Float64}}()
    @test !isdiscrete(m5)
    @test iscontinuous(m5)
    @test isnormalized(m5)
    @test domaintype(m5) == SVector{2,Float64}
    @test support(m5) == FullSpace{SVector{2,Float64}}()
    @test weight(m5, SVector(0.0,0.0)) ≈ 1/(2pi)
end
