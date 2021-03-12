
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

default_atol() = default_atol(Float64)
default_rtol() = default_rtol(Float64)

default_atol(::Type{T}) where {T} = zero(T)
default_rtol(::Type{T}) where {T} = sqrt(eps(T))
default_rtol(x) = default_rtol(typeof(x))

default_maxevals() = 10^4

"An adaptive quadrature strategy"
struct QuadAdaptive{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int
end

QuadAdaptive(; atol = default_atol(), rtol = default_rtol(atol), maxevals = default_maxevals()) = QuadAdaptive(atol, rtol, maxevals)
QuadAdaptive(atol, rtol, maxevals) = QuadAdaptive(promote(atol,rtol)..., maxevals)
QuadAdaptive(atol::T, rtol::T, maxevals) where {T} = QuadAdaptive{T}(atol, rtol, maxevals)

QuadAdaptive{T}(atol = default_atol(T), rtol = default_rtol(atol)) where {T} =
    QuadAdaptive{T}(atol, rtol, default_maxevals())

"Adaptive quadrature using quadgk"
struct Q_quadgk{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int
end

Q_quadgk(; atol = default_atol(), rtol = default_rtol(atol), maxevals = default_maxevals()) = Q_quadgk(atol, rtol, maxevals)
Q_quadgk(atol, rtol, maxevals) = Q_quadgk(promote(atol,rtol)..., maxevals)
Q_quadgk(atol::T, rtol::T, maxevals) where {T} = Q_quadgk{T}(atol, rtol, maxevals)

Q_quadgk{T}(atol = default_atol(T), rtol = default_rtol(atol)) where {T} =
    Q_quadgk{T}(atol, rtol, default_maxevals())

Q_quadgk(qs::QuadAdaptive{T}) where {T} = Q_quadgk{T}(qs.atol, qs.rtol, qs.maxevals)


"Adaptive quadrature using hcubature"
struct Q_hcubature{T} <: AdaptiveStrategy{T}
    atol    ::  T
    rtol    ::  T
    maxevals    ::  Int
end

Q_hcubature(; atol = default_atol(), rtol = default_rtol(atol), maxevals = default_maxevals()) =
    Q_hcubature(atol, rtol, maxevals)
Q_hcubature(atol, rtol, maxevals) = Q_hcubature(promote(atol,rtol)..., maxevals)
Q_hcubature(atol::T, rtol::T, maxevals) where {T} = Q_hcubature{T}(atol, rtol, maxevals)

Q_hcubature{T}(atol = default_atol(T), rtol = default_rtol(atol)) where {T} =
    Q_hcubature{T}(atol, rtol, default_maxevals())

Q_hcubature(qs::QuadAdaptive{T}) where {T} = Q_hcubature{T}(qs.atol, qs.rtol, qs.maxevals)



"Use a quadrature rule adapted to the measure."
struct BestRule <: QuadratureStrategy
    n   ::  Int
end


"A quadrature rule associated with a domain."
abstract type DomainRule{T} <: QuadratureStrategy end

points(q::DomainRule) = q.x
weights(q::DomainRule) = q.w

abstract type IntervalRule{T} <: DomainRule{T} end

"Apply a quadrature rule defined on `[-1,1]`."
struct ChebyshevIntervalRule{T} <: IntervalRule{T}
    x       ::  Vector{T}
    w       ::  Vector{T}
end
domain(q::ChebyshevIntervalRule{T}) where {T} = ChebyshevInterval{T}()


"Apply a fixed quadrature rule defined on `[0,1]`."
struct UnitIntervalRule{T} <: IntervalRule{T}
    x       ::  Vector{T}
    w       ::  Vector{T}
end
domain(q::UnitIntervalRule{T}) where {T} = UnitInterval{T}()


"Apply a fixed quadrature rule defined on the halfline `[0,∞)`."
struct HalfLineRule{T} <: DomainRule{T}
    x       ::  Vector{T}
    w       ::  Vector{T}
end
domain(q::HalfLineRule{T}) where {T} = HalfLine{T}()

"Apply a fixed quadrature rule defined on the real line."
struct RealLineRule{T} <: DomainRule{T}
    x       ::  Vector{T}
    w       ::  Vector{T}
end
domain(q::RealLineRule{T}) where {T} = FullSpace{T}()


struct GenericDomainRule{T} <: DomainRule{T}
    domain  ::  Domain{T}
    x       ::  Vector{T}
    w       ::  Vector{T}
end
domain(q::GenericDomainRule) = q.domain


Q_GaussLegendre(args...) = ChebyshevIntervalRule(gausslegendre(args...)...)
Q_GaussJacobi(args...) = ChebyshevIntervalRule(gaussjacobi(args...)...)
Q_GaussLaguerre(args...) = HalfLineRule(gausslaguerre(args...)...)
Q_GaussHermite(args...) = RealLineRule(gausshermite(args...)...)


"""
Turn the interval rule into a composite quadrature rule with `L` segments on
the same interval.
"""
function composite_rule(L::Int, rule::ChebyshevIntervalRule{T}) where {T}
    x = points(rule)
    w = weights(rule)
    grid = range(-one(T), one(T), length=L+1)
    h = step(grid)
    composite_x = vcat([grid[i] .+ (x .+ 1)*h/2 for i in 1:L]...)
    composite_w = vcat([w*h/2 for i in 1:L]...)
    ChebyshevIntervalRule(composite_x, composite_w)
end

composite_rule(L::Int, n::Int, ::Type{T} = Float64) where {T} =
    composite_rule(L, Q_GaussLegendre(T, n))


"""
Construct a graded quadrature rule towards the left of the unit interval.

The rule is a composite quadrature rule on intervals of the form
`σ^(k+1):σ^k`. The number of points on interval `i` is `⌈n*i*μ⌉`,
where `μ` is an optional parameter.
"""
function graded_rule_left(σ, n, μ = 1, δ = 100eps(σ), nmax = 30)
    k = floor(Int, log(δ)/log(σ))
    graded_mesh = σ.^(k:-1:0)
    T = typeof(σ)
    x = T[]
    w = T[]
    for i in 1:k
        q = Q_GaussLegendre(T, min(ceil(Int, n*i*μ),nmax))
        ai = graded_mesh[i]
        bi = graded_mesh[i+1]
        h = bi-ai
        x = vcat(x, ai .+ (points(q) .+ 1)*h/2)
        w = vcat(w, weights(q)*h/2)
    end
    UnitIntervalRule(x, w)
end

"""
Construct a graded quadrature rule towards the right of the unit interval.

See: `graded_rule_left`.
"""
function graded_rule_right(σ, n, args...)
    Q = graded_rule_left(σ, n, args...)
    x = points(Q)
    w = weights(Q)
    xrev = reverse(1 .- x)
    wrev = reverse(w)
    UnitIntervalRule(xrev, wrev)
end
