
export integral,
    integrate


"""
The suggested quadrature strategy based on the `integral` arguments.

By default, adaptive quadrature is chosen.
"""
suggestedstrategy(args...) = QuadAdaptive()

returntype(integrand, ::Type{S}) where {S} = Base.Core.Compiler.return_type(integrand, (S,))

zero_error(x) = zero_error(typeof(x))
zero_error(::Type{T}) where {T} = zero(prectype(T))
unknown_error(T) = -one(prectype(T))
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
tolerance(::Type{T}) where {T} = tolerance(prectype(T))


"""
Compute an integral of the given integrand on the given domain,
or with the given measure.

Example:
```
integral(cos, 0.0..1.0)
```
"""
function integral end

"Like integral, but also returns an error estimate (if applicable)."
function integrate end

function integral(args...)
    I, E = integrate(args...)
    I
end

# associate a domain with a measure: use its support
_domain(μ::Measure) = support(μ)

# associate a measure with a domain
# we try to avoid memory allocations, hence LebesgeeSpace{T} is the default
_measure(domain::Domain{T}) where {T} = Lebesgue{T}()
# these cases have a known allocation-free lebesgue measure
_measure(domain::ChebyshevInterval) = lebesguemeasure(domain)
_measure(domain::UnitInterval) = lebesguemeasure(domain)
# sometimes the T of an interval is an integer (e.g. in 0..1)
_measure(domain::AbstractInterval{T}) where {T} = Lebesgue{float(T)}()


# Process the arguments until there is a domain, a measure and a property object.
process_arguments(measure::Measure) =
    process_arguments(_domain(measure), measure)
process_arguments(domain::Domain, measure::Measure) =
    process_arguments(domain, measure, NoProperty())
process_arguments(measure::Measure, domain::Domain) =
    process_arguments(domain, measure, NoProperty())
process_arguments(measure::Measure, prop::Property) =
    process_arguments(_domain(measure), measure, prop)
process_arguments(domain::Domain{T}) where {T} =
    process_arguments(domain, _measure(domain))
process_arguments(domain::Domain{T}, prop::Property) where {T} =
    process_arguments(domain, _measure(domain), prop)
# all good now
process_arguments(domain::Domain, measure::Measure, prop::Property) =
    (domain, measure, prop)

process_arguments(args...) = error("Arguments to integral or integrate functions not understood.")

integrate(qs::ChebyshevIntervalRule{T}, integrand) where {T} =
    integrate(qs, integrand, ChebyshevInterval{T}())
integrate(qs::HalfLineRule{T}, integrand) where {T} =
    integrate(qs, integrand, HalfLine{T}())
integrate(qs::RealLineRule{T}, integrand) where {T} =
    integrate(qs, integrand, FullSpace{T}())

## The use of generators
integrate(gen::Base.Generator, args...) =
    integrate(process_generator(gen)..., args...)
# - catch something like f(x) for x in domain
process_generator(gen::Base.Generator{<:Domain}) = (gen.f, gen.iter)
# - catch something like f(x,y) for x in domain1, y in domain2
process_generator(gen::Base.Generator{<:Base.Iterators.ProductIterator}) =
    process_generator(gen, gen.iter.iterators)
function process_generator(gen, iterators::Tuple{Vararg{Domain}})
    domain = productdomain(iterators...)
    dims = map(dimension, iterators)
    @show dims
    @show typeof(iterators)
    (gen.f, domain)
end

integrate(integrand, args...) =
    integrate(integrand, process_arguments(args...)...)

integrate(qs::QuadratureStrategy, integrand, args...) =
    integrate(qs, integrand, process_arguments(args...)...)

integrate(integrand, domain::Domain, measure::Measure, prop::Property) =
    integrate(suggestedstrategy(domain, measure, prop), integrand, domain, measure, prop)


# The process is as follows:
# - The argument list is completed (missing domain or measure are added)
#   -> integrate_prop is invoked
# - integrate_prop processes properties using dispatch
#   -> integrate_measure is invoked
# - integrate_measure processes measures using dispatch
#   -> integrate_domain is invoked
# - integrate_domain processes domains using dispatch
#   -> apply_quad is invoked
# - If apply_quad is not intercepted, a fallback routine is invoked.



struct Integrand{F}
    fun ::  F
end

(∘)(f::Integrand, g::Function) = integrand_compose(f.fun, g)
integrand_compose(f, g) = Integrand(t -> f(g(t)))
(∘)(f::Integrand, g::IdentityMap) = f

(*)(f::IdentityMap, g::Integrand) = g
(*)(f, g::Integrand) = integrand_times(f, g.fun)
integrand_times(f::Function, g) = Integrand(t -> f(t)*g(t))
integrand_times(c::Number, g) = Integrand(t -> c*g(t))

# go to step P once we have the correct signature
integrate(qs::QuadratureStrategy, integrand::Integrand, domain::Domain, measure::Measure, prop::Property) =
    integrate_prop(qs, integrand, domain, measure, prop)
integrate(qs::QuadratureStrategy, integrand, domain::Domain, measure::Measure, prop::Property) =
    integrate_prop(qs, Integrand(integrand), domain, measure, prop)

# STEP P: dispatch on the property
integrate_prop(qs, integrand, domain, measure, prop) =
    integrate_measure(qs, integrand, domain, measure, prop)

# STEP M: dispatch on the measure
integrate_measure(qs, integrand, domain, measure, prop) =
    integrate_domain(qs, integrand, domain, measure, prop)

# STEP D: dispatch on the domain
integrate_domain(qs, integrand, domain, measure, prop) =
    integrate_done(qs, integrand, domain, measure, prop)

integrate_done(qs, integrand::Integrand, domain, measure, prop) =
    integrate_done(qs, integrand.fun, domain, measure, prop)
integrate_done(qs, integrand::Function, domain, measure, prop) =
    select_quad(qs, integrand, domain, measure, prop)

# Selection: use a suitable quadrature strategy

# Generic adaptive quadrature in 1D invokes QuadGK, everywhere else it uses hcubature
select_quad(qs::QuadAdaptive, integrand, domain::Domain{T}, measure, prop) where {T<:Number} =
    apply_quad(Q_quadgk(qs), integrand, domain, measure, prop)
select_quad(qs::QuadAdaptive, integrand, domain::Domain{T}, measure, prop) where {T} =
    apply_quad(Q_hcubature(qs), integrand, domain, measure, prop)
select_quad(qs, integrand, domain, measure, prop) =
    apply_quad(qs, integrand, domain, measure, prop)

# FINAL STEP: invoke apply_quad, with a fallback if necessary
apply_quad(qs, integrand, domain, measure, prop) =
    fallback_integrate(qs, integrand, domain, measure, prop)

apply_quad(qs, integrand, domain::ProductDomain, measure, prop) =
    apply_productquad(qs, integrand, domain, measure, prop, components(domain)...)

apply_productquad(qs, integrand, domain, measure, prop, domains...) =
    fallback_integrate(qs, integrand, domain, measure, prop)

fallback_integrate(qs, integrand, domain, measure, prop) =
    error("Don't know how to integrate on $domain with strategy $(qs).")


# For quadgk, we only know how to compute intervals
function apply_quad(qs::Q_quadgk, integrand, domain::AbstractInterval, measure::LebesgueMeasure, prop)
    if isempty(domain)
        zero_result(integrand, prectype(domain))
    else
        quadgk(integrand, extrema(domain)...; atol = qs.atol, rtol = qs.rtol, maxevals = qs.maxevals)
    end
end

function apply_productquad(qs::Q_hcubature, integrand, domain, measure, prop, domains::AbstractInterval...)
    if islebesguemeasure(measure)
        apply_hcubature(qs,integrand, domains...)
    else
        error("HCubature is invoked, but not with a Lebesgue measure: $(measure)")
    end
end

function apply_hcubature(qs::Q_hcubature, integrand, domains::AbstractInterval...)
    a = map(leftendpoint, domains)
    b = map(rightendpoint, domains)
    hcubature(integrand, a, b; rtol=qs.rtol,atol=qs.atol, maxevals = qs.maxevals)
end

# Given a rule on [-1,1] and a different interval, we map the rule
function apply_quad(qs::ChebyshevIntervalRule, integrand, domain::AbstractInterval, measure::LebesgueMeasure, prop)
    a = leftendpoint(domain)
    b = rightendpoint(domain)
    D = (b-a)/2
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(a + (x[i]+1)*D)*D for i in 1:length(x))
    z, unknown_error(z)
end

function apply_quad(qs::UnitIntervalRule, integrand, domain::AbstractInterval, measure::LebesgueMeasure, prop)
    a = leftendpoint(domain)
    b = rightendpoint(domain)
    D = b-a
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(a + x[i]*D)*D for i in 1:length(x))
    z, unknown_error(z)
end

function apply_quad(qs::IntervalRule, integrand, interval::AbstractInterval, measure::LebesgueMeasure, prop)
    a, b = extrema(interval)
    A, B = extrema(domain(qs))
    m = mapto(A..B, a..b)
    J = (b-a)/(B-A)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(m(x[i]))*J for i in 1:length(x))
    z, unknown_error(z)
end

# Make sure that a half line rule is applied to a half line integral
function apply_quad(qs::HalfLineRule, integrand, domain::HalfLine, measure::LebesgueMeasure, prop)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end

# Make sure that a real line rule is applied to a real line integral
function apply_quad(qs::RealLineRule, integrand, domain::FullSpace, measure::LebesgueMeasure, prop)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end

# Any other quadrature rule, we just apply it
function apply_quad(qs::GenericDomainRule, integrand, domain, measure, prop)
    @assert domain(qs) == domain
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end
