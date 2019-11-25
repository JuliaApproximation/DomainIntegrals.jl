
abstract type QuadratureStrategy end

abstract type AdaptiveStrategy <: QuadratureStrategy end

abstract type FixedRuleStrategy <: QuadratureStrategy end

abstract type RuleFamilyStrategy <: QuadratureStrategy end

struct QuadAdaptive <: AdaptiveStrategy
end

struct Q_quadgk <: AdaptiveStrategy
end

struct Q_hcubature <: AdaptiveStrategy
end


struct Q_GaussLegendre{T} <: FixedRuleStrategy
    x   ::  Array{T,1}
    w   ::  Array{T,1}
end

Q_GaussLegendre(n::Int = 4) = Q_GaussLegendre{Float64}(n)
function Q_GaussLegendre{Float64}(n::Int)
    x, w = FastGaussQuadrature.gausslegendre(n)
    Q_GaussLegendre(x, w)
end

function Q_GaussLegendre{T}(n::Int) where {T}
    x, w = GaussQuadrature.legendre(T, n)
    Q_GaussLegendre{T}(x, w)
end
