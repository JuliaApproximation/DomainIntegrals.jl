## Basic weight functions

export LebesgueMeasure,
    Lebesgue,
    LebesgueDomain,
    lebesguemeasure,
    LegendreWeight,
    JacobiWeight,
    ChebyshevWeight,
    ChebyshevTWeight,
    ChebyshevUWeight,
    LaguerreWeight,
    HermiteWeight,
    GaussianWeight,
    DiracWeight,
    point


####################
# Lebesgue measures
####################


"Supertype of Lebesgue measures"
abstract type LebesgueMeasure{T} <: Weight{T} end

unsafe_weightfun(μ::LebesgueMeasure, x) = one(codomaintype(μ))

islebesguemeasure(m::Measure) = false
islebesguemeasure(m::LebesgueMeasure) = true

"The Lebesgue measure on the space `FullSpace{T}`."
struct Lebesgue{T} <: LebesgueMeasure{T}
end

Lebesgue() = Lebesgue{Float64}()
similar(μ::Lebesgue, ::Type{T}) where {T} = Lebesgue{T}()


"The Lebesgue measure on the unit interval `[0,1]`."
struct LebesgueUnit{T} <: LebesgueMeasure{T}
end

LebesgueUnit() = LebesgueUnit{Float64}()
similar(μ::LebesgueUnit, ::Type{T}) where {T <: Real} = Lebesgue{T}()
support(μ::LebesgueUnit{T}) where {T} = UnitInterval{T}()

isnormalized(μ::LebesgueUnit) = true

"Lebesgue measure supported on a general domain."
struct LebesgueDomain{T} <: LebesgueMeasure{T}
    domain  ::  Domain{T}
end

similar(μ::LebesgueDomain, ::Type{T}) where {T} = LebesgueDomain{T}(μ.domain)
support(m::LebesgueDomain) = m.domain



###############
# Dirac weight
###############

"A continuous Dirac measure at a point"
struct DiracWeight{T} <: Weight{T}
    point   ::  T
end

similar(μ::DiracWeight, ::Type{T}) where {T} = DiracWeight{T}(μ.point)

point(μ::DiracWeight) = μ.point
support(μ::DiracWeight) = Point(μ.point)
isnormalized(μ::DiracWeight) = true
unsafe_weightfun(μ::DiracWeight, x) = convert(codomaintype(μ), Inf)



###################################
# Classical orthogonal polynomials
###################################


"The Legendre weight is a Lebesgue measure on the interval `[-1,1]`."
struct LegendreWeight{T} <: LebesgueMeasure{T}
end
LegendreWeight() = LegendreWeight{Float64}()

similar(μ::LegendreWeight, ::Type{T}) where {T <: Real} = LegendreWeight{T}()
support(μ::LegendreWeight{T}) where {T} = ChebyshevInterval{T}()


"The Jacobi weight on the interval `[-1,1]`."
struct JacobiWeight{T} <: Weight{T}
    α   ::  T
    β   ::  T

    JacobiWeight{T}(α = zero(T), β = zero(T)) where {T} = new(α, β)
end
JacobiWeight() = JacobiWeight{Float64}()
JacobiWeight(α, β) = JacobiWeight(promote(α, β)...)
JacobiWeight(α::T, β::T) where {T<:AbstractFloat} = JacobiWeight{T}(α, β)
JacobiWeight(α::N, β::N) where {N<:Number} = JacobiWeight(float(α), float(β))

jacobi_weightfun(x, α, β) = (1+x)^α * (1-x)^β

similar(μ::JacobiWeight, ::Type{T}) where {T <: Real} = JacobiWeight{T}(μ.α, μ.β)
support(μ::JacobiWeight{T}) where {T} = ChebyshevInterval{T}()
unsafe_weightfun(μ::JacobiWeight, x) = jacobi_weightfun(x, μ.α, μ.β)

jacobi_α(μ::JacobiWeight) = μ.α
jacobi_β(μ::JacobiWeight) = μ.β


"""
The `Chebyshev` or `ChebyshevT` weight is the measure on `[-1,1]` with the
Chebyshev weight function `w(x) = 1/√(1-x^2)`.
"""
struct ChebyshevTWeight{T} <: Weight{T}
end
ChebyshevTWeight() = ChebyshevTWeight{Float64}()

const ChebyshevWeight = ChebyshevTWeight

chebyshev_weight_firstkind(x) = 1/sqrt(1-x^2)

similar(μ::ChebyshevTWeight, ::Type{T}) where {T <: Real} = ChebyshevTWeight{T}()
support(μ::ChebyshevTWeight{T}) where {T} = ChebyshevInterval{T}()
unsafe_weightfun(μ::ChebyshevTWeight, x) = chebyshev_weight_firstkind(x)

weightfunction(μ::ChebyshevTWeight) = chebyshev_weight_firstkind


"""
The ChebyshevU weight is the measure on `[-1,1]` with the Chebyshev weight
function of the second kind `w(x) = √(1-x^2).`
"""
struct ChebyshevUWeight{T} <: Weight{T}
end
ChebyshevUWeight() = ChebyshevUWeight{Float64}()

chebyshev_weight_secondkind(x) = sqrt(1-x^2)

similar(μ::ChebyshevUWeight, ::Type{T}) where {T <: Real} = ChebyshevUWeight{T}()
support(μ::ChebyshevUWeight{T}) where {T} = ChebyshevInterval{T}()
unsafe_weightfun(μ::ChebyshevUWeight, x) = chebyshev_weight_secondkind(x)

weightfunction(μ::ChebyshevUWeight) = chebyshev_weight_secondkind


convert(::Type{JacobiWeight}, μ::LegendreWeight{T}) where {T} =
    JacobiWeight{T}(0, 0)
convert(::Type{JacobiWeight}, μ::ChebyshevTWeight{T}) where {T} =
    JacobiWeight{T}(-one(T)/2, -one(T)/2)
convert(::Type{JacobiWeight}, μ::ChebyshevUWeight{T}) where {T} =
    JacobiWeight{T}(one(T)/2, one(T)/2)

function convert(::Type{ChebyshevTWeight}, μ::JacobiWeight{T}) where {T}
    (jacobi_α(μ) ≈ -one(T)/2 && jacobi_β(μ) ≈ -one(T)/2) || throw(InexactError(:convert, ChebyshevTWeight, μ))
    ChebyshevTWeight{T}()
end

function convert(::Type{ChebyshevUWeight}, μ::JacobiWeight{T}) where {T}
    (jacobi_α(μ) ≈ one(T)/2 && jacobi_β(μ) ≈ one(T)/2) || throw(InexactError(:convert, ChebyshevUWeight, μ))
    ChebyshevUWeight{T}()
end

function convert(::Type{LegendreWeight}, μ::JacobiWeight{T}) where {T}
    (μ.α ≈ 0 && μ.β ≈ 0) || throw(InexactError(:convert, LegendreWeight, μ))
    LegendreWeight{T}()
end

jacobi_α(μ::LegendreWeight{T}) where {T} = zero(T)
jacobi_β(μ::LegendreWeight{T}) where {T} = zero(T)
jacobi_α(μ::ChebyshevTWeight{T}) where {T} = -one(T)/2
jacobi_β(μ::ChebyshevTWeight{T}) where {T} = -one(T)/2
jacobi_α(μ::ChebyshevUWeight{T}) where {T} = one(T)/2
jacobi_β(μ::ChebyshevUWeight{T}) where {T} = one(T)/2



"The generalised Laguerre measure on the halfline `[0,∞)`."
struct LaguerreWeight{T} <: Weight{T}
    α   ::  T

    LaguerreWeight{T}(α = zero(T)) where {T} = new(α)
end
LaguerreWeight() = LaguerreWeight{Float64}()
LaguerreWeight(α::T) where {T<:AbstractFloat} = LaguerreWeight{T}(α)
LaguerreWeight(α) = LaguerreWeight(float(α))

laguerre_weightfun(x, α) = exp(-x) * x^α

similar(μ::LaguerreWeight, ::Type{T}) where {T <: Real} = LaguerreWeight{T}(μ.α)
support(μ::LaguerreWeight{T}) where {T} = HalfLine{T}()
isnormalized(m::LaguerreWeight) = m.α == 0

unsafe_weightfun(μ::LaguerreWeight, x) = laguerre_weightfun(x, μ.α)

laguerre_α(μ::LaguerreWeight) = μ.α


"The Hermite measure with weight exp(-x^2) on the real line."
struct HermiteWeight{T} <: Weight{T}
end
HermiteWeight() = HermiteWeight{Float64}()

hermite_weightfun(x) = exp(-x^2)

similar(μ::HermiteWeight, ::Type{T}) where {T <: Real} = HermiteWeight{T}()
unsafe_weightfun(μ::HermiteWeight, x) = hermite_weightfun(x)

weightfunction(μ::HermiteWeight) = hermite_weightfun


"The Gaussian measure with weight exp(-|x|^2/2)."
struct GaussianWeight{T} <: Weight{T}
end
GaussianWeight() = GaussianWeight{Float64}()

gaussian_weightfun(x, ::Type{T} = prectype(x)) where {T} =
    1/(2*convert(T, pi))^(length(x)/2) * exp(-norm(x)^2)

similar(μ::GaussianWeight, ::Type{T}) where {T} = GaussianWeight{T}()
isnormalized(μ::GaussianWeight) = true
unsafe_weightfun(μ::GaussianWeight, x) = gaussian_weightfun(x, prectype(μ))



"The Lebesgue measure associated with the given domain"
lebesguemeasure(domain::UnitInterval{T}) where {T} = LebesgueUnit{T}()
lebesguemeasure(domain::ChebyshevInterval{T}) where {T} = LegendreWeight{T}()
lebesguemeasure(domain::FullSpace{T}) where {T} = Lebesgue{T}()
lebesguemeasure(domain::Domain{T}) where {T} = LebesgueDomain{T}(domain)
