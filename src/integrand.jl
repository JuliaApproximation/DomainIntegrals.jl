
"Supertype of an integrand which evaluates to `T`."
abstract type Integrand{T} end

zero_result(f) = zero_result(f, Float64)
zero_result(integrand::Integrand{T}) where {T} = zero_result(integrand, T)

zero_result(f, ::Type{T}) where T = zero(T), zero(real(T))

convert(::Type{Integrand{T}}, f::Integrand{T}) where {T} = f
convert(::Type{Integrand{T}}, f::Integrand{S}) where {S,T} =
    similar_integrand(f, T)

Base.complex(f::Integrand{T}) where {T} =
    convert(Integrand{complex(T)}, F)

(∘)(f::Integrand, g) = integrand_compose(f, g, f.fun)
(∘)(f::Integrand, g::IdentityMap) = f

(*)(f, g::Integrand{T}) where {T} = integrand_times(f, g, g.fun)
(*)(f::Map, g::Integrand{T}) where {T} = integrand_times(t->applymap(f,t), g, g.fun)
(*)(f::ConstantMap, g::Integrand) = mapconstant(f) * g


"Representation of an integrand function."
struct FunIntegrand{T,F} <: Integrand{T}
    fun ::  F
end

FunIntegrand{T}(fun::F) where {T,F} = FunIntegrand{T,F}(fun)
FunIntegrand{T}(F::FunIntegrand) where {T} = FunIntegrand{T}(F.fun)

similar_integrand(F::FunIntegrand, ::Type{T}) where {T} = FunIntegrand{T}(F)

integrand_compose(f::FunIntegrand{T}, g, fun) where {T} = FunIntegrand{T}(t -> fun(g(t)))

integrand_times(f::Function, g::FunIntegrand{T}, fun) where {T} =
    FunIntegrand{T}(t -> f(t)*fun(t))
integrand_times(c::S, g::FunIntegrand{T}, fun) where {T,S<:Number} =
    FunIntegrand{promote_type(S,T)}(t -> c*fun(t))

(f::FunIntegrand)(x...) = f.fun(x...)


"An integrand that logs any transformations which are applied to it."
struct TransformationLogIntegrand{T} <: Integrand{T}
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
