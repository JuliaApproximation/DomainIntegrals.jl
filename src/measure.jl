
"""
A measure on a set (or `Domain{T}`) is a function that assigns a non-negative
number to subsets of that set. The support of a `Measure{T}` is a `Domain{T}`.

We use measures mostly to define weighted integrals, i.e., we consider
measures `μ` for which a weight function `w(x)` exists satisfying
`dμ = w(x) dx`.
"""
abstract type Measure{T} end

export support
"Return the support of the measure?"
function support end

export domaintype
domaintype(μ::Measure) = domaintype(typeof(μ))
domaintype(::Type{<:Measure{T}}) where {T} = T

export codomaintype
"What is the codomain type of the measure?"
codomaintype(μ::Measure) = codomaintype(typeof(μ))
codomaintype(::Type{<:Measure{T}}) where {T} = prectype(T)

export weight
"The weight function associated with the measure."
weight(μ::Measure{T}, x) where {T} = weight(μ, convert(T, x))
weight(μ::Measure{T}, x::T) where {T} =
    x ∈ support(μ) ? unsafe_weight(μ, x) : zero(codomaintype(μ))

"Supertype of Lebesgue measures"
abstract type AbstractLebesgueMeasure{T} <: Measure{T} end

unsafe_weight(μ::AbstractLebesgueMeasure, x) = one(codomaintype(μ))

export LebesgueMeasure
"The Lebesgue measure on the space `FullSpace{T}`."
struct LebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

LebesgueMeasure() = LebesgueMeasure{Float64}()

support(μ::LebesgueMeasure{T}) where {T} = FullSpace{T}()

export UnitLebesgueMeasure
"The Lebesgue measure on the unit interval `[0,1]`."
struct UnitLebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

UnitLebesgueMeasure() = UnitLebesgueMeasure{Float64}()

support(μ::UnitLebesgueMeasure{T}) where {T} = UnitInterval{T}()

export LegendreMeasure
"The Lebesgue measure on the interval `[-1,1]`."
struct LegendreMeasure{T} <: AbstractLebesgueMeasure{T}
end
LegendreMeasure() = LegendreMeasure{Float64}()

support(μ::LegendreMeasure{T}) where {T} = ChebyshevInterval{T}()


export JacobiMeasure
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


export LaguerreMeasure
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

export HermiteMeasure
"The Hermite measure with weight exp(-x^2) on the real line."
struct HermiteMeasure{T} <: Measure{T}
end
HermiteMeasure() = HermiteMeasure{Float64}()

support(μ::HermiteMeasure{T}) where {T} = FullSpace{T}()

unsafe_weight(μ::HermiteMeasure, x) = exp(-x^2)


export DiracMeasure
"A point measure"
struct DiracMeasure{T} <: Measure{T}
    point   ::  T
end

export point
point(μ::DiracMeasure) = μ.point
support(μ::DiracMeasure) = Point(μ.point)

unsafe_weight(μ::DiracMeasure, x) = convert(codomaintype(μ), Inf)
