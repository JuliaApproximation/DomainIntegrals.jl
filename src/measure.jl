
"""
A measure on a set (or `Domain{T}`) is a function that assigns a non-negative
number to subsets of that set. The support of a `Measure{T}` is a `Domain{T}`.

We use measures mostly to define weighted integrals, i.e., we consider
measures `μ` for which a weight function `w(x)` exists satisfying
`dμ = w(x) dx`.
"""
abstract type Measure{T} end

domaintype(μ::Measure) = domaintype(typeof(μ))
domaintype(::Type{<:Measure{T}}) where {T} = T

codomaintype(μ::Measure) = codomaintype(typeof(μ))
codomaintype(::Type{<:Measure{T}}) where {T} = prectype(T)

weight(μ::Measure{T}, x) where {T} = weight(μ, convert(T, x))
weight(μ::Measure{T}, x::T) where {T} =
    x ∈ support(μ) ? unsafe_weight(μ, x) : zero(rangetype(μ))


abstract type AbstractLebesgueMeasure{T} <: Measure{T} end

unsafe_weight(μ::AbstractLebesgueMeasure, x) = one(codomaintype(μ))

"The Lebesgue measure on the space `FullSpace{T}`."
struct LebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

support(μ::LebesgueMeasure{T}) where {T} = FullSpace{T}()

"The Lebesgue measure on the unit interval `[0,1]`."
struct UnitLebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

support(μ::UnitLebesgueMeasure{T}) where {T} = UnitInterval{T}()

"The Lebesgue measure on the interval `[-1,1]`."
struct LegendreMeasure{T} <: AbstractLebesgueMeasure{T}
end

support(μ::LegendreMeasure{T}) where {T} = ChebyshevInterval{T}()


struct DiracMeasure{T} <: Measure{T}
    point   ::  T
end

point(μ::DiracMeasure) = μ.point
support(μ::DiracMeasure) = Point(μ.point)

unsafe_weight(μ::DiracMeasure, x) = convert(codomaintype(μ), Inf)
