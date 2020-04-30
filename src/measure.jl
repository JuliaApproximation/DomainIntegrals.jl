
export support,
    domaintype,
    codomaintype,
    isdiscrete,
    iscontinuous,
    isnormalized,
    weight,
    weightfunction,
    points,
    weights,
    discrete_weight,
    AbstractLebesgueMeasure,
    LebesgueMeasure,
    UnitLebesgueMeasure,
    DomainLebesgueMeasure,
    LegendreMeasure,
    JacobiMeasure,
    ChebyshevTMeasure,
    ChebyshevMeasure,
    ChebyshevUMeasure,
    LaguerreMeasure,
    HermiteMeasure,
    GaussianMeasure,
    DiracMeasure,
    point,
    lebesguemeasure


"Supertype of measures. The support of an `AbstractMeasure{T}` is a `Domain{T}`."
abstract type AbstractMeasure{T} end

domaintype(μ::AbstractMeasure) = domaintype(typeof(μ))
domaintype(::Type{<:AbstractMeasure{T}}) where {T} = T

"What is the codomain type of the measure?"
codomaintype(μ::AbstractMeasure) = codomaintype(typeof(μ))
codomaintype(::Type{<:AbstractMeasure{T}}) where {T} = prectype(T)

prectype(::Type{<:AbstractMeasure{T}}) where {T} = prectype(T)
numtype(::Type{<:AbstractMeasure{T}}) where {T} = numtype(T)

"Is the measure normalized?"
isnormalized(μ::AbstractMeasure) = false

convert(::Type{AbstractMeasure{T}}, μ::AbstractMeasure{T}) where {T} = μ
convert(::Type{AbstractMeasure{T}}, μ::AbstractMeasure{S}) where {S,T} = similar(μ, T)


"""
A `Measure{T}` is a continuous measure that is defined in terms of a
weightfunction: `dμ = w(x) dx`.
"""
abstract type Measure{T} <: AbstractMeasure{T} end

"""
A `DiscreteMeasure` is defined in terms of a discrete set of points and an
associated weight vector.
"""
abstract type DiscreteMeasure{T} <: AbstractMeasure{T} end

"Is the measure discrete?"
isdiscrete(μ::Measure) = false
isdiscrete(μ::DiscreteMeasure) = true

"Is the measure continuous?"
iscontinuous(μ::Measure) = true
iscontinuous(μ::DiscreteMeasure) = false


## Support for continuous measures

"Return the support of the measure"
support(μ::Measure{T}) where {T} = FullSpace{T}()

"Evaluate the weight function associated with the measure."
function weight(μ::AbstractMeasure{T}, x::S) where {S,T}
    U = promote_type(S,T)
    weight(convert(AbstractMeasure{U}, μ), convert(U, x))
end
weight(μ::AbstractMeasure{T}, x::T) where {T} = weight1(μ, x)

weight1(μ::AbstractMeasure, x) =
    x ∈ support(μ) ? unsafe_weight(μ, x) : zero(codomaintype(μ))

weightfunction(m::AbstractMeasure) = x->weight(m, x)
unsafe_weightfunction(m::AbstractMeasure) = x->unsafe_weight(m, x)

## Support for discrete measures

points(μ::DiscreteMeasure) = μ.points
weights(μ::DiscreteMeasure) = μ.weights

support(μ::DiscreteMeasure) = points(μ)

function discrete_weight(μ::DiscreteMeasure, i)
    # Perform a bounds check and invoke unsafe_discrete_weight,
    # so that concrete measures may implement e.g. an on-the-fly formula for
    # the weights without bounds checking
    @boundscheck checkbounds(μ, i)
    unsafe_discrete_weight(μ, i)
end
checkbounds(μ::DiscreteMeasure, i) = checkbounds(points(μ), i)

function unsafe_discrete_weight(μ::DiscreteMeasure, i)
    @inbounds weights(μ)[i]
end

isnormalized(μ::DiscreteMeasure) = sum(weights(μ)) ≈ 1

length(μ::DiscreteMeasure) = length(points(μ))


## Basic measures

"Supertype of Lebesgue measures"
abstract type AbstractLebesgueMeasure{T} <: Measure{T} end

unsafe_weight(μ::AbstractLebesgueMeasure, x) = one(codomaintype(μ))

"The Lebesgue measure on the space `FullSpace{T}`."
struct LebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

LebesgueMeasure() = LebesgueMeasure{Float64}()
similar(μ::LebesgueMeasure, ::Type{T}) where {T} = LebesgueMeasure{T}()


"The Lebesgue measure on the unit interval `[0,1]`."
struct UnitLebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
end

UnitLebesgueMeasure() = UnitLebesgueMeasure{Float64}()
similar(μ::UnitLebesgueMeasure, ::Type{T}) where {T <: Real} = LebesgueMeasure{T}()
support(μ::UnitLebesgueMeasure{T}) where {T} = UnitInterval{T}()

isnormalized(μ::UnitLebesgueMeasure) = true

"Lebesgue measure supported on a general domain."
struct DomainLebesgueMeasure{T} <: AbstractLebesgueMeasure{T}
    domain  ::  Domain{T}
end

similar(μ::DomainLebesgueMeasure, ::Type{T}) where {T} = DomainLebesgueMeasure{T}(μ.domain)
support(m::DomainLebesgueMeasure) = m.domain



#################
# Applications
#################

"A point measure"
struct DiracMeasure{T} <: DiscreteMeasure{T}
    point   ::  T
end

similar(μ::DiracMeasure, ::Type{T}) where {T} = DiracMeasure{T}(μ.point)

point(μ::DiracMeasure) = μ.point
support(μ::DiracMeasure) = Point(μ.point)
isnormalized(μ::DiracMeasure) = true
unsafe_weight(μ::DiracMeasure, x) = convert(codomaintype(μ), Inf)



## Some widely used measures associated with orthogonal polynomials follow


"The Lebesgue measure on the interval `[-1,1]`."
struct LegendreMeasure{T} <: AbstractLebesgueMeasure{T}
end
LegendreMeasure() = LegendreMeasure{Float64}()

similar(μ::LegendreMeasure, ::Type{T}) where {T <: Real} = LegendreMeasure{T}()
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

similar(μ::JacobiMeasure, ::Type{T}) where {T <: Real} = JacobiMeasure{T}(μ.α, μ.β)
support(μ::JacobiMeasure{T}) where {T} = ChebyshevInterval{T}()
unsafe_weight(μ::JacobiMeasure, x) = (1+x)^μ.α * (1-x)^μ.β


"""
The `Chebyshev` or `ChebyshevT` measure is the measure on `[-1,1]` with the
Chebyshev weight `w(x) = 1/√(1-x^2)`.
"""
struct ChebyshevTMeasure{T} <: Measure{T}
end
ChebyshevTMeasure() = ChebyshevTMeasure{Float64}()

const ChebyshevMeasure = ChebyshevTMeasure

similar(μ::ChebyshevTMeasure, ::Type{T}) where {T <: Real} = ChebyshevTMeasure{T}()
support(μ::ChebyshevTMeasure{T}) where {T} = ChebyshevInterval{T}()
unsafe_weight(μ::ChebyshevTMeasure, x) = 1/sqrt(1-x^2)


"""
The ChebyshevU measure is the measure on `[-1,1]` with the Chebyshev weight
of the second kind `w(x) = √(1-x^2).`
"""
struct ChebyshevUMeasure{T} <: Measure{T}
end
ChebyshevUMeasure() = ChebyshevUMeasure{Float64}()

similar(μ::ChebyshevUMeasure, ::Type{T}) where {T <: Real} = ChebyshevUMeasure{T}()
support(μ::ChebyshevUMeasure{T}) where {T} = ChebyshevInterval{T}()
unsafe_weight(μ::ChebyshevUMeasure, x) = sqrt(1-x^2)


convert(::Type{JacobiMeasure}, μ::LegendreMeasure{T}) where {T} =
    JacobiMeasure{T}(0, 0)
convert(::Type{JacobiMeasure}, μ::ChebyshevTMeasure{T}) where {T} =
    JacobiMeasure{T}(-one(T)/2, -one(T)/2)
convert(::Type{JacobiMeasure}, μ::ChebyshevUMeasure{T}) where {T} =
    JacobiMeasure{T}(one(T)/2, one(T)/2)

function convert(::Type{ChebyshevTMeasure}, μ::JacobiMeasure{T}) where {T}
    (μ.α ≈ -one(T)/2 && μ.β ≈ -one(T)/2) || throw(InexactError(:convert, ChebyshevTMeasure, μ))
    ChebyshevTMeasure{T}()
end

function convert(::Type{ChebyshevUMeasure}, μ::JacobiMeasure{T}) where {T}
    (μ.α ≈ one(T)/2 && μ.β ≈ one(T)/2) || throw(InexactError(:convert, ChebyshevUMeasure, μ))
    ChebyshevUMeasure{T}()
end

function convert(::Type{LegendreMeasure}, μ::JacobiMeasure{T}) where {T}
    (μ.α ≈ 0 && μ.β ≈ 0) || throw(InexactError(:convert, LegendreMeasure, μ))
    LegendreMeasure{T}()
end


"The generalised Laguerre measure on the halfline `[0,∞)`."
struct LaguerreMeasure{T} <: Measure{T}
    α   ::  T

    LaguerreMeasure{T}(α = zero(T)) where {T} = new(α)
end
LaguerreMeasure() = LaguerreMeasure{Float64}()
LaguerreMeasure(α::T) where {T<:AbstractFloat} = LaguerreMeasure{T}(α)
LaguerreMeasure(α) = LaguerreMeasure(float(α))

similar(μ::LaguerreMeasure, ::Type{T}) where {T <: Real} = LaguerreMeasure{T}(μ.α)
support(μ::LaguerreMeasure{T}) where {T} = HalfLine{T}()
isnormalized(m::LaguerreMeasure) = m.α == 0
unsafe_weight(μ::LaguerreMeasure, x) = exp(-x) * x^μ.α


"The Hermite measure with weight exp(-x^2) on the real line."
struct HermiteMeasure{T} <: Measure{T}
end
HermiteMeasure() = HermiteMeasure{Float64}()

similar(μ::HermiteMeasure, ::Type{T}) where {T <: Real} = HermiteMeasure{T}()
unsafe_weight(μ::HermiteMeasure, x) = exp(-x^2)


"The Gaussian measure with weight exp(-|x|^2/2)."
struct GaussianMeasure{T} <: Measure{T}
end
GaussianMeasure() = GaussianMeasure{Float64}()

similar(μ::GaussianMeasure, ::Type{T}) where {T} = GaussianMeasure{T}()
isnormalized(μ::GaussianMeasure) = true
unsafe_weight(μ::GaussianMeasure, x) = 1/(2*convert(prectype(μ), pi))^(length(x)/2) * exp(-norm(x)^2)



"The Lebesgue measure associated with the given domain"
lebesguemeasure(domain::UnitInterval{T}) where {T} = UnitLebesgueMeasure{T}()
lebesguemeasure(domain::ChebyshevInterval{T}) where {T} = LegendreMeasure{T}()
lebesguemeasure(domain::FullSpace{T}) where {T} = LebesgueMeasure{T}()
lebesguemeasure(domain::Domain{T}) where {T} = DomainLebesgueMeasure{T}(domain)
