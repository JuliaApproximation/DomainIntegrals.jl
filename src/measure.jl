
export support,
    domaintype,
    codomaintype,
    isdiscrete,
    iscontinuous,
    isnormalized,
    weight,
    LebesgueMeasure,
    UnitLebesgueMeasure,
    LegendreMeasure,
    JacobiMeasure,
    LaguerreMeasure,
    HermiteMeasure,
    DiracMeasure,
    point

"""
Supertype of all measures.

The support of an `AbstractMeasure{T}` is a `Domain{T}`.
"""
abstract type AbstractMeasure{T} end

"""
A `Measure{T}` is a continuous measure that is defined in terms of a
weightfunction: `dμ = w(x) dx`.
"""
abstract type Measure{T} <: AbstractMeasure{T} end

"""
A `DiscreteMeasure` is defined in terms of a discrete grid and an
associated weight vector.
"""
abstract type DiscreteMeasure{T} <: Measure{T} end

"Return the support of the measure"
support(μ::Measure{T}) where {T} = FullSpace{T}()

"Is the measure discrete?"
isdiscrete(μ::Measure) = false
isdiscrete(μ::DiscreteMeasure) = true

"Is the measure continuous?"
iscontinuous(μ::Measure) = true
iscontinuous(μ::DiscreteMeasure) = false


"Is the measure normalized?"
isnormalized(μ::AbstractMeasure) = false

domaintype(μ::AbstractMeasure) = domaintype(typeof(μ))
domaintype(::Type{<:AbstractMeasure{T}}) where {T} = T

"What is the codomain type of the measure?"
codomaintype(μ::AbstractMeasure) = codomaintype(typeof(μ))
codomaintype(::Type{<:AbstractMeasure{T}}) where {T} = prectype(T)

"Evaluate the weight function associated with the measure."
weight(μ::Measure{T}, x) where {T} = weight(μ, convert(T, x))
weight(μ::Measure{T}, x::T) where {T} =
    x ∈ support(μ) ? unsafe_weight(μ, x) : zero(codomaintype(μ))


#################
# Basic measures
#################

"Supertype of Lebesgue measures"
abstract type AbstractLebesgueMeasure{T} <: Measure{T} end

unsafe_weight(μ::AbstractLebesgueMeasure, x) = one(codomaintype(μ))

"The Lebesgue measure on the space `FullSpace{T}`."
struct LebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

LebesgueMeasure() = LebesgueMeasure{Float64}()


"The Lebesgue measure on the unit interval `[0,1]`."
struct UnitLebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

UnitLebesgueMeasure() = UnitLebesgueMeasure{Float64}()

support(μ::UnitLebesgueMeasure{T}) where {T} = UnitInterval{T}()


#################
# Applications
#################


"A point measure"
struct DiracMeasure{T} <: DiscreteMeasure{T}
    point   ::  T
end

point(μ::DiracMeasure) = μ.point
support(μ::DiracMeasure) = Point(μ.point)

isnormalized(μ::DiracMeasure) = true

unsafe_weight(μ::DiracMeasure, x) = convert(codomaintype(μ), Inf)



## Some widely used measures associated with orthogonal polynomials follow


"The Lebesgue measure on the interval `[-1,1]`."
struct LegendreMeasure{T} <: AbstractLebesgueMeasure{T}
end
LegendreMeasure() = LegendreMeasure{Float64}()

support(μ::LegendreMeasure{T}) where {T} = ChebyshevInterval{T}()


"The Jacobi measure on the interval `[-1,1]`."
struct JacobiMeasure{T} <: Measure{T}
    α   ::  T
    β   ::  T

    JacobiMeasure{T}(α = zero(T), β = zero(T)) where {T} = new(α, β)
end
JacobiMeasure() = JacobiMeasure{Float64}()
JacobiMeasure(α, β) = JacobiMeasure(promote(α, β)...)
JacobiMeasure(α::T, β::T) where {T<:AbstractFloat} = JacobiMeasure{T}(α, β)
JacobiMeasure(α::N, β::N) where {N<:Number} = JacobiMeasure(float(α), float(β))

support(μ::JacobiMeasure{T}) where {T} = ChebyshevInterval{T}()

unsafe_weight(μ::JacobiMeasure, x) = (1+x)^μ.α * (1-x)^μ.β


"The generalised Laguerre measure on the halfline `[0,∞)`."
struct LaguerreMeasure{T} <: Measure{T}
    α   ::  T

    LaguerreMeasure{T}(α = zero(T)) where {T} = new(α)
end
LaguerreMeasure() = LaguerreMeasure{Float64}()
LaguerreMeasure(α::T) where {T<:AbstractFloat} = LaguerreMeasure{T}(α)
LaguerreMeasure(α) = LaguerreMeasure(float(α))

support(μ::LaguerreMeasure{T}) where {T} = HalfLine{T}()

unsafe_weight(μ::LaguerreMeasure, x) = exp(-x) * x^μ.α


"The Hermite measure with weight exp(-x^2) on the real line."
struct HermiteMeasure{T} <: Measure{T}
end
HermiteMeasure() = HermiteMeasure{Float64}()

support(μ::HermiteMeasure{T}) where {T} = FullSpace{T}()

unsafe_weight(μ::HermiteMeasure, x) = exp(-x^2)
