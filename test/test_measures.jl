
function test_measures()
    m1 = Lebesgue()
    @test !isdiscrete(m1)
    @test iscontinuous(m1)
    @test !isnormalized(m1)
    @test domaintype(m1) == Float64
    @test DomainIntegrals.unsafe_weightfun(m1, 0.4) == 1

    m2 = DomainIntegrals.LebesgueUnit()
    @test !isdiscrete(m2)
    @test iscontinuous(m2)
    @test isnormalized(m2)
    @test domaintype(m2) == Float64
    @test support(m2) == UnitInterval()
    @test weightfun(m2, 0.4) == 1
    @test weightfun(m2, 1.4) == 0
    @test m2 == dx(UnitInterval())

    m3 = LegendreWeight()
    @test !isdiscrete(m3)
    @test iscontinuous(m3)
    @test !isnormalized(m3)
    @test domaintype(m3) == Float64
    @test support(m3) == ChebyshevInterval()
    @test weightfun(m3, 0.4) == 1
    @test weightfun(m3, -0.4) == 1
    @test weightfun(m3, 1.4) == 0
    @test weightfun(m3, -1.4) == 0
    @test jacobi_α(m3) == 0
    @test jacobi_β(m3) == 0

    x = 0.5
    m4 = DiracWeight(x)
    @test !isdiscrete(m4)
    @test iscontinuous(m4)
    @test isnormalized(m4)
    @test domaintype(m4) == Float64
    @test support(m4) == Point(x)
    @test point(m4) == x
    @test weightfun(m4, x) == Inf
    @test weightfun(m4, x+1) == 0
    @test weightfun(m4, big(x)) == Inf
    @test weightfun(m4, big(x+1)) == 0

    m5 = GaussianWeight{SVector{2,Float64}}()
    @test !isdiscrete(m5)
    @test iscontinuous(m5)
    @test isnormalized(m5)
    @test domaintype(m5) == SVector{2,Float64}
    @test support(m5) == FullSpace{SVector{2,Float64}}()
    @test weightfun(m5, SVector(0.0,0.0)) ≈ 1/(2pi)
    @test weightfun(m5, SVector(0,0)) ≈ 1/(2pi)
    @test weightfun(m5, SVector(big(0),big(0))) ≈ 1/(2pi)
    @test weightfun(m5, SVector(big(0.0),big(0.0))) ≈ 1/(2pi)

    m6 = LaguerreWeight(0.0)
    @test !isdiscrete(m6)
    @test iscontinuous(m6)
    @test isnormalized(m6)
    @test !isnormalized(LaguerreWeight(0.1))
    @test domaintype(m6) == Float64
    @test support(m6) == HalfLine()
    @test weightfun(m6, 0.4) == exp(-0.4)
    @test weightfun(m6, -0.4) == 0
    @test weightfun(m6, big(0.4)) == exp(-big(0.4))

    m7 = ChebyshevTWeight()
    @test support(m7) isa ChebyshevInterval
    @test !isdiscrete(m7)
    @test iscontinuous(m7)
    @test !isnormalized(m7)
    @test domaintype(m7) == Float64
    @test jacobi_α(m7) == -1/2
    @test jacobi_β(m7) == -1/2
    @test m7 == JacobiWeight(-1/2, -1/2)

    m8 = ChebyshevUWeight()
    @test support(m8) isa ChebyshevInterval
    @test !isdiscrete(m8)
    @test iscontinuous(m8)
    @test !isnormalized(m8)
    @test domaintype(m8) == Float64
    @test jacobi_α(m8) == 1/2
    @test jacobi_β(m8) == 1/2
    @test m8 == JacobiWeight(1/2, 1/2)

    m9 = JacobiWeight(0.4, 0.7)
    @test support(m9) == (-1..1)
    @test !isdiscrete(m9)
    @test iscontinuous(m9)
    @test !isnormalized(m9)
    @test domaintype(m8) == Float64
    @test jacobi_α(m9) == 0.4
    @test jacobi_β(m9) == 0.7
    @test weightfun(m9, 0.22) ≈ (1-0.22)^0.4 * (1+0.22)^0.7
    @test_throws AssertionError JacobiWeight(-2, 0)
    @test_throws AssertionError JacobiWeight(0, -1)

    m10 = UltrasphericalWeight(0.3)
    @test support(m10) == (-1..1)
    @test !isdiscrete(m10)
    @test iscontinuous(m10)
    @test !isnormalized(m10)
    @test domaintype(m10) == Float64
    @test jacobi_α(m10) == 0.3-0.5
    @test jacobi_β(m10) == 0.3-0.5
    @test weightfun(m10, 0.22) ≈ (1-0.22^2)^(0.3-0.5)
end

function test_discrete_measures()
    x = [0.5, 1.0]
    w = [0.2, 0.8]
    μ = DomainIntegrals.GenericDiscreteWeight(x, w)
    @test isnormalized(μ)
    @test !isuniform(μ)
    x2, w2 = μ
    @test norm(x-x2) ≈ 0
    @test norm(w-w2) ≈ 0
end
