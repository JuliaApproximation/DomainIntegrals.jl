
"Supertype of an integrand which evaluates to `T`."
abstract type AbstractIntegrand{T} end

zero_result(f) = 0.0, 0.0
zero_result(integrand::AbstractIntegrand{T}) where {T} = zero(T), zero(real(T))

convert(::Type{AbstractIntegrand{T}}, f::AbstractIntegrand{T}) where {T} = f
convert(::Type{AbstractIntegrand{T}}, f::AbstractIntegrand{S}) where {S,T} =
    similar_integrand(f, T)

Base.complex(f::AbstractIntegrand{T}) where {T} =
    convert(AbstractIntegrand{complex(T)}, F)

(∘)(f::AbstractIntegrand, g) = integrand_compose(f, g, f.fun)
(∘)(f::AbstractIntegrand, g::IdentityMap) = f

(*)(f, g::AbstractIntegrand{T}) where {T} = integrand_times(f, g, g.fun)
(*)(f::Map, g::AbstractIntegrand{T}) where {T} = integrand_times(t->applymap(f,t), g, g.fun)
(*)(f::ConstantMap, g::AbstractIntegrand) = constant(f) * g


"Representation of an integrand function."
struct Integrand{T,F} <: AbstractIntegrand{T}
    fun ::  F
end

Integrand{T}(fun::F) where {T,F} = Integrand{T,F}(fun)
Integrand{T}(F::Integrand) where {T} = Integrand{T}(F.fun)

similar_integrand(F::Integrand, ::Type{T}) where {T} = Integrand{T}(F)

integrand_compose(f::Integrand{T}, g, fun) where {T} = Integrand{T}(t -> fun(g(t)))

integrand_times(f::Function, g::Integrand{T}, fun) where {T} =
    Integrand{T}(t -> f(t)*fun(t))
integrand_times(c::S, g::Integrand{T}, fun) where {T,S<:Number} =
    Integrand{promote_type(S,T)}(t -> c*fun(t))

(f::Integrand)(x...) = f.fun(x...)


"An integrand that logs any transformations which are applied to it."
struct TransformationLogIntegrand{T} <: AbstractIntegrand{T}
    fun
    prefactor
end

TransformationLogIntegrand{T}() where {T} = TransformationLogIntegrand{T}(IdentityMap{T}())
TransformationLogIntegrand{T}(fun) where {T} = TransformationLogIntegrand{T}(fun, UnityMap{T}())

similar_integrand(F::TransformationLogIntegrand, ::Type{T}) where {T} =
    TransformationLogIntegrand{T}(F)

integrand_compose(f::TransformationLogIntegrand{T}, g, fun::IdentityMap) where {T} =
    TransformationLogIntegrand{T}(g, f.prefactor)

integrand_compose(f::TransformationLogIntegrand{T}, g, fun) where {T} =
    TransformationLogIntegrand{T}(fun ∘ g, f.prefactor)

log_times(f::Function, g::UnityMap) = f
log_times(f::Function, g) = t -> f(t)*g(t)
log_times(c::Number, g::UnityMap{T}) where {T} = ConstantMap{T}(c)

integrand_times(f, g::TransformationLogIntegrand{T}, fun) where {T} =
    TransformationLogIntegrand{T}(fun, log_times(f, g.prefactor))


"Wrap a function and log any point at which it is evaluated."
struct LoggingFunction{T,F}
    fun     ::  F
    points  ::  Vector{T}
end

LoggingFunction{T}(f::F) where {T,F} = LoggingFunction{T,F}(f, Vector{T}())

function (f::LoggingFunction)(x)
    push!(f.points, x)
    f.fun(x)
end
