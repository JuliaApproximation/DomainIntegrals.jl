
export QuadAdaptive,
    BestRule,
    Q_GaussLegendre,
    Q_GaussJacobi,
    Q_GaussLaguerre,
    Q_GaussHermite


abstract type QuadratureStrategy end

abstract type AdaptiveStrategy{T} <: QuadratureStrategy end

prectype(::Type{<:AdaptiveStrategy{T}}) where {T} = prectype(T)


"An adaptive quadrature strategy"
struct QuadAdaptive{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int

    QuadAdaptive{T}(atol = zero(T), rtol = sqrt(eps(T)), maxevals = 10^4) where {T} = new(atol, rtol, maxevals)
end

QuadAdaptive() = QuadAdaptive{Float64}()

"Adaptive quadrature using quadgk"
struct Q_quadgk{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int

    Q_quadgk{T}(atol = zero(T), rtol = sqrt(eps(T)), maxevals = 10^4) where {T} = new(atol, rtol, maxevals)
end

Q_quadgk() = Q_quadgk{Float64}()
Q_quadgk(qs::QuadAdaptive{T}) where {T} = Q_quadgk{T}(qs.atol, qs.rtol, qs.maxevals)


"Adaptive quadrature using hcubature"
struct Q_hcubature{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int

    Q_hcubature{T}(atol = zero(T), rtol = sqrt(eps(T)), maxevals = 10^4) where {T} = new(atol, rtol, maxevals)
end

Q_hcubature() = Q_hcubature{Float64}()
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
