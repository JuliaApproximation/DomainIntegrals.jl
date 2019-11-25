
"""
The suggested quadrature strategy based on the arguments supplied to
the `integral` function.

By default, adaptive quadrature is chosen.
"""
suggestedstrategy(args...) = QuadAdaptive()

returntype(integrand, ::Type{S}) where {S} = Base.Core.Compiler.return_type(integrand, (S,))

errortype(x) = errortype(typeof(x))
errortype(::Type{T}) where {T <: AbstractFloat} = T
errortype(::Type{Complex{T}}) where {T} = errortype(T)

unknown_error(T) = -one(errortype(T))
zero_error(T) = zero(errortype(T))

function zero_result(integrand, ::Type{S}) where {S}
    T = returntype(integrand, S)
    zero(T), zero_error(T)
end

tolerance(d::Domain{T}) where {T} = tolerance(T)
tolerance(::Type{T}) where {T<:AbstractFloat} = 10eps(T)
tolerance(::Type{Complex{T}}) where {T} = tolerance(T)
tolerance(::Type{SVector{N,T}}) where {N,T} = tolerance(T)

"""
Compute an integral of the given integrand on the given domain,
or with the given measure.

Example:
```
integral(cos, 0.0..1.0)
```
"""
integral(integrand) = error("Please specify an integration domain or an integration measure.")
function integral(integrand, arg1, args...)
    I, E = quadrature(suggestedstrategy(arg1, args...), integrand, arg1, args...)
    I
end
function integral(qs::QuadratureStrategy, integrand, arg1, args...)
    I, E = quadrature(qs, integrand, arg1, args...)
    I
end



# The process is as follows:
# - The argument list is completed (missing domain or measure are added)
#   -> quadrature_s is invoked
# - quadrature_s processes singularities using dispatch
#   -> quadrature_m is invoked
# - quadrature_m processes measures using dispatch
#   -> quadrature_d is invoked
# - quadrature_d processes domains using dispatch
#   -> apply_quad is invoked
# - If apply_quad is not intercepted, a fallback routine is invoked.

# Completion of interface:
# - the default domain is the support of the measure (if given)
# - the default measure for Domain{T} is LebesgueMeasure{T}
# - the default singularity is no singularity
quadrature(qs::QuadratureStrategy, integrand, measure::Measure) =
    quadrature(qs, integrand, support(measure), measure)
quadrature(qs::QuadratureStrategy, integrand, measure::Measure, sing::Singularity) =
    quadrature(qs, integrand, support(measure), measure, sing)
quadrature(qs::QuadratureStrategy, integrand, domain::Domain{T}) where {T} =
    quadrature(qs, integrand, domain, LebesgueMeasure{T}())
quadrature(qs::QuadratureStrategy, integrand, domain::Domain{T}, sing::Singularity) where {T} =
    quadrature(qs, integrand, domain, LebesgueMeasure{T}(), sing)
quadrature(qs::QuadratureStrategy, integrand, domain::Domain, measure::Measure) =
    quadrature(qs, integrand, domain, measure, NoSingularity())
quadrature(integrand, arg1, args...) =
    quadrature(suggestedstrategy(arg1, args...), integrand, arg1, args...)

# go to step S once we have the correct signature
quadrature(qs::QuadratureStrategy, integrand, domain::Domain, measure::Measure, sing::Singularity) =
    quadrature_s(qs, integrand, domain, measure, sing)

# STEP S: dispatch on the singularity
quadrature_s(qs, integrand, domain, measure, sing) =
    quadrature_m(qs, integrand, domain, measure, sing)

# STEP M: dispatch on the measure
quadrature_m(qs, integrand, domain, measure, sing) =
    quadrature_d(qs, integrand, domain, measure, sing)

# STEP D: dispatch on the domain
quadrature_d(qs, integrand, domain, measure, sing) =
    select_quad(qs, integrand, domain, measure, sing)


# Selection: use a suitable quadrature strategy

# Generic adaptive quadrature in 1D invokes QuadGK, everywhere else it uses hcubature
select_quad(::QuadAdaptive, integrand, domain::Domain{T}, measure, sing) where {T<:Number} =
    apply_quad(Q_quadgk(), integrand, domain, measure, sing)
select_quad(::QuadAdaptive, integrand, domain::Domain{T}, measure, sing) where {T} =
    apply_quad(Q_hcubature(), integrand, domain, measure, sing)
select_quad(qs, integrand, domain, measure, sing) =
    apply_quad(qs, integrand, domain, measure, sing)

# FINAL STEP: invoke apply_quad, with a fallback if necessary
apply_quad(qs, integrand, domain, measure, sing) =
    fallback_quadrature(qs, integrand, domain, measure, sing)

apply_quad(qs, integrand, domain::ProductDomain, measure, sing) =
    apply_productquad(qs, integrand, domain, measure, sing, elements(domain)...)

apply_productquad(qs, integrand, domain, measure, sing, domains...) =
    fallback_quadrature(qs, integrand, domain, measure, sing)

fallback_quadrature(qs, integrand, domain, measure, sing) =
    error("Don't know how to integrate on $domain.")


# For quadgk, we only know how to compute intervals
apply_quad(::Q_quadgk, integrand, domain::AbstractInterval, measure::AbstractLebesgueMeasure, sing) =
    quadgk(integrand, extrema(domain)...)

# hcubature works for rectangles with the Lebesgue measure
apply_productquad(::Q_hcubature, integrand, domain, measure::AbstractLebesgueMeasure, sing, domains::AbstractInterval...) =
    apply_hcubature(integrand, domains...)


function apply_hcubature(integrand, domains::AbstractInterval...)
    a = map(leftendpoint, domains)
    b = map(rightendpoint, domains)
    hcubature(integrand, a, b)
end

function apply_quad(qs::Q_GaussLegendre, integrand, domain::AbstractInterval, measure::AbstractLebesgueMeasure, sing)
    a = leftendpoint(domain)
    b = rightendpoint(domain)
    D = (b-a)/2
    T = typeof(integrand(a))
    z = zero(T)
    n = length(qs.x)
    for i in 1:n
        z += sum(qs.w[i]*integrand(a + (qs.x[i]+1)*D)*D)
    end
    z, unknown_error(z)
end
