
export integral,
    quadrature

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

zero_error(T) = zero(errortype(T))
unknown_error(T) = -one(errortype(T))
unknown_error(z::Number) = -one(z)

function zero_result(integrand, ::Type{S}) where {S}
    T = returntype(integrand, S)
    if T == Any
        zero(S), zero_error(S)
    else
        zero(T), zero_error(T)
    end
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
function integral(qs::QuadratureStrategy, integrand, args...)
    I, E = quadrature(qs, integrand, args...)
    I
end

"Like integral, but also returns an error estimate (if applicable)."
function quadrature end


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
quadrature(qs::QuadratureStrategy, integrand, measure::AbstractMeasure) =
    quadrature(qs, integrand, support(measure), measure)
quadrature(qs::QuadratureStrategy, integrand, measure::AbstractMeasure, sing::Singularity) =
    quadrature(qs, integrand, support(measure), measure, sing)
quadrature(qs::QuadratureStrategy, integrand, domain::Domain{T}) where {T} =
    quadrature(qs, integrand, domain, LebesgueMeasure{T}())
quadrature(qs::QuadratureStrategy, integrand, domain::Domain{T}, sing::Singularity) where {T} =
    quadrature(qs, integrand, domain, LebesgueMeasure{T}(), sing)
quadrature(qs::QuadratureStrategy, integrand, measure::AbstractMeasure, domain::Domain) =
    quadrature(qs, integrand, domain, measure)
quadrature(qs::QuadratureStrategy, integrand, domain::Domain, measure::AbstractMeasure) =
    quadrature(qs, integrand, domain, measure, NoSingularity())
quadrature(qs::ChebyshevIntervalRule{T}, integrand) where {T} =
    quadrature(qs, integrand, ChebyshevInterval{T}())
quadrature(qs::HalfLineRule{T}, integrand) where {T} =
    quadrature(qs, integrand, HalfLine{T}())
quadrature(qs::RealLineRule{T}, integrand) where {T} =
    quadrature(qs, integrand, FullSpace{T}())

quadrature(::QuadratureStrategy, arg1, args...) =
    error("Incorrect type or number of arguments?")
quadrature(integrand, arg1, args...) =
    quadrature(suggestedstrategy(arg1, args...), integrand, arg1, args...)

# go to step S once we have the correct signature
quadrature(qs::QuadratureStrategy, integrand, domain::Domain, measure::AbstractMeasure, sing::Singularity) =
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
select_quad(qs::QuadAdaptive, integrand, domain::Domain{T}, measure, sing) where {T<:Number} =
    apply_quad(Q_quadgk(qs), integrand, domain, measure, sing)
select_quad(qs::QuadAdaptive, integrand, domain::Domain{T}, measure, sing) where {T} =
    apply_quad(Q_hcubature(qs), integrand, domain, measure, sing)
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
    error("Don't know how to integrate on $domain with strategy $(qs).")


# For quadgk, we only know how to compute intervals
function apply_quad(qs::Q_quadgk, integrand, domain::AbstractInterval, measure::AbstractLebesgueMeasure, sing)
    if isempty(domain)
        zero_result(integrand, prectype(qs))
    else
        quadgk(integrand, extrema(domain)...; atol = qs.atol, rtol = qs.rtol, maxevals = qs.maxevals)
    end
end

function apply_productquad(qs::Q_hcubature, integrand, domain, measure, sing, domains::AbstractInterval...)
    if islebesguemeasure(measure)
        apply_hcubature(qs,integrand, domains...)
    else
        error("HCubature is invoked, but not with a Lebesgue measure: $(measure)")
    end
end

function apply_hcubature(qs::Q_hcubature,integrand, domains::AbstractInterval...)
    a = map(leftendpoint, domains)
    b = map(rightendpoint, domains)
    hcubature(integrand, a, b; rtol=qs.rtol,atol=qs.atol, maxevals = qs.maxevals)
end

# Given a rule on [-1,1] and a different interval, we map the rule
function apply_quad(qs::ChebyshevIntervalRule, integrand, domain::AbstractInterval, measure::AbstractLebesgueMeasure, sing)
    a = leftendpoint(domain)
    b = rightendpoint(domain)
    D = (b-a)/2
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(a + (x[i]+1)*D)*D for i in 1:length(x))
    z, unknown_error(z)
end

function apply_quad(qs::UnitIntervalRule, integrand, domain::AbstractInterval, measure::AbstractLebesgueMeasure, sing)
    a = leftendpoint(domain)
    b = rightendpoint(domain)
    D = b-a
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(a + x[i]*D)*D for i in 1:length(x))
    z, unknown_error(z)
end

function apply_quad(qs::IntervalRule, integrand, interval::AbstractInterval, measure::AbstractLebesgueMeasure, sing)
    a, b = extrema(interval)
    A, B = extrema(domain(qs))
    m = interval_map(A, B, a, b)
    J = (b-a)/(B-A)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(m(x[i]))*J for i in 1:length(x))
    z, unknown_error(z)
end

# Make sure that a half line rule is applied to a half line integral
function apply_quad(qs::HalfLineRule, integrand, domain::HalfLine, measure::AbstractLebesgueMeasure, sing)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end

# Make sure that a real line rule is applied to a real line integral
function apply_quad(qs::RealLineRule, integrand, domain::FullSpace, measure::AbstractLebesgueMeasure, sing)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end

# Any other quadrature rule, we just apply it
function apply_quad(qs::GenericDomainRule, integrand, domain, measure, sing)
    @assert domain(qs) == domain
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end
