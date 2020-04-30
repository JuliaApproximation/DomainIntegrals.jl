
export QuadAdaptive,
    Q_quadgk,
    Q_hcubature,
    BestRule,
    Q_GaussLegendre,
    Q_GaussJacobi,
    Q_GaussLaguerre,
    Q_GaussHermite


abstract type QuadratureStrategy end

abstract type AdaptiveStrategy{T} <: QuadratureStrategy end

prectype(::Type{<:AdaptiveStrategy{T}}) where {T} = prectype(T)

def_atol() = def_atol(Float64)
def_rtol() = def_rtol(Float64)

def_atol(::Type{T}) where {T} = zero(T)
def_rtol(::Type{T}) where {T} = sqrt(eps(T))
def_rtol(x) = def_rtol(typeof(x))

def_maxevals() = 10^4

"An adaptive quadrature strategy"
struct QuadAdaptive{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int
end

QuadAdaptive(; atol = def_atol(), rtol = def_rtol(atol), maxevals = def_maxevals()) = QuadAdaptive(atol, rtol, maxevals)
QuadAdaptive(atol, rtol, maxevals) = QuadAdaptive(promote(atol,rtol)..., maxevals)
QuadAdaptive(atol::T, rtol::T, maxevals) where {T} = QuadAdaptive{T}(atol, rtol, maxevals)

QuadAdaptive{T}(atol = def_atol(T), rtol = def_rtol(atol)) where {T} =
    QuadAdaptive{T}(atol, rtol, def_maxevals())

"Adaptive quadrature using quadgk"
struct Q_quadgk{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int
end

Q_quadgk(; atol = def_atol(), rtol = def_rtol(atol), maxevals = def_maxevals()) = Q_quadgk(atol, rtol, maxevals)
Q_quadgk(atol, rtol, maxevals) = Q_quadgk(promote(atol,rtol)..., maxevals)
Q_quadgk(atol::T, rtol::T, maxevals) where {T} = Q_quadgk{T}(atol, rtol, maxevals)

Q_quadgk{T}(atol = def_atol(T), rtol = def_rtol(atol)) where {T} =
    Q_quadgk{T}(atol, rtol, def_maxevals())

Q_quadgk(qs::QuadAdaptive{T}) where {T} = Q_quadgk{T}(qs.atol, qs.rtol, qs.maxevals)


"Adaptive quadrature using hcubature"
struct Q_hcubature{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int
end

Q_hcubature(; atol = 0.0, rtol = sqrt(eps(atol)), maxevals = 10^4) = Q_hcubature(atol, rtol, maxevals)
Q_hcubature(atol, rtol, maxevals) = Q_hcubature(promote(atol,rtol)..., maxevals)
Q_hcubature(atol::T, rtol::T, maxevals) where {T} = Q_hcubature{T}(atol, rtol, maxevals)

Q_hcubature{T}(atol = def_atol(T), rtol = def_rtol(atol)) where {T} =
    Q_hcubature{T}(atol, rtol, def_maxevals())

Q_hcubature(qs::QuadAdaptive{T}) where {T} = Q_hcubature{T}(qs.atol, qs.rtol, qs.maxevals)


"Apply a fixed quadrature rule"
abstract type FixedRule <: QuadratureStrategy end

points(s::FixedRule) = s.x
weights(s::FixedRule) = s.w

"Use a quadrature rule adapted to the measure."
struct BestRule <: QuadratureStrategy
    n   ::  Int
end


"Apply a fixed quadrature rule defined on `[-1,1]`."
struct FixedRuleInterval{T} <: FixedRule
    x       ::  Vector{T}
    w       ::  Vector{T}
end

"Apply a fixed quadrature rule defined on the halfline `[0,∞)`."
struct FixedRuleHalfLine{T} <: FixedRule
    x       ::  Vector{T}
    w       ::  Vector{T}
end

"Apply a fixed quadrature rule defined on the real line."
struct FixedRuleRealLine{T} <: FixedRule
    x       ::  Vector{T}
    w       ::  Vector{T}
end

Q_GaussLegendre(args...) = FixedRuleInterval(gausslegendre(args...)...)
Q_GaussJacobi(args...) = FixedRuleInterval(gaussjacobi(args...)...)
Q_GaussLaguerre(args...) = FixedRuleHalfLine(gausslaguerre(args...)...)
Q_GaussHermite(args...) = FixedRuleRealLine(gausshermite(args...)...)
