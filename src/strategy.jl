
abstract type QuadratureStrategy end

abstract type AdaptiveStrategy <: QuadratureStrategy end


export QuadAdaptive
"An adaptive quadrature strategy"
struct QuadAdaptive <: AdaptiveStrategy
end

"Adaptive quadrature using quadgk"
struct Q_quadgk <: AdaptiveStrategy
end

"Adaptive quadrature using hcubature"
struct Q_hcubature <: AdaptiveStrategy
end

"Apply a fixed quadrature rule"
abstract type FixedRule <: QuadratureStrategy end

points(s::FixedRule) = s.x
weights(s::FixedRule) = s.w

export BestRule
"Use a quadrature rule adapted to the measure."
struct BestRule <: QuadratureStrategy
    n   ::  Int
end

export Q_GaussLegendre,
    Q_GaussJacobi,
    Q_GaussLaguerre,
    Q_GaussHermite

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
