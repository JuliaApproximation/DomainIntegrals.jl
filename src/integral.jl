
export integral,
    integrate


"""
The suggested quadrature strategy based on the `integral` arguments.

By default, adaptive quadrature is chosen.
"""
suggestedstrategy(domain, args...) = QuadAdaptive{prectype(DomainRef(domain))}()

returntype(integrand, ::Type{S}) where {S} = Base.Core.Compiler.return_type(integrand, (S,))

zero_error(x) = zero_error(typeof(x))
zero_error(::Type{T}) where {T} = zero(prectype(T))
unknown_error(T) = -one(prectype(T))
unknown_error(z::Number) = -one(z)

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
associated_domain(μ::Measure) = support(μ)
associated_domain(μ::DiscreteWeight{T}) where {T} = DummyDiscreteDomain{T}()

# associate a measure with a domain
associated_measure(domain) = Lebesgue{domaineltype(domain)}()

promote_domain_and_measure(domain::Domain{T}, measure::Measure{T}) where {T} =
    domain, measure
function promote_domain_and_measure(domain::Domain{S}, measure::Measure{T}) where {S,T}
    U = promote_type(S,T)
    convert(Domain{U}, domain), convert(Measure{U}, measure)
end


# Process the arguments until there is a domain (not necessarily of Domain type),
# a measure and zero or more property objects.
process_arguments(measure::Measure, properties::Property...) =
    process_arguments(associated_domain(measure), measure, properties...)
process_arguments(domain, properties::Property...) =
    process_arguments(checkdomain(domain), associated_measure(domain), properties...)

function process_arguments(domain::Domain{S}, measure::Measure{T}, properties::Property...) where {S,T}
    d, m = promote_domain_and_measure(domain, measure)
    process_arguments(d, m, properties...)
end
process_arguments(domain::Domain{T}, measure::Measure{T}, properties::Property...) where {T} =
    (domain, measure, properties...)
process_arguments(domain, measure::Measure{T}, properties::Property...) where {T} =
    (domain, measure, properties...)

process_arguments(args...) = error("Arguments to integral or integrate functions not understood.")

integrate(qs::ChebyshevIntervalRule{T}, integrand) where {T} =
    integrate(qs, integrand, ChebyshevInterval{T}())
integrate(qs::HalfLineRule{T}, integrand) where {T} =
    integrate(qs, integrand, HalfLine{T}())
integrate(qs::RealLineRule{T}, integrand) where {T} =
    integrate(qs, integrand, RealLine{T}())

## The use of generators
integrate(gen::Base.Generator, args...) =
    integrate(process_generator(gen)..., args...)
# - catch something like f(x) for x in domain
process_generator(gen::Base.Generator{<:AnyDomain}) = (gen.f, gen.iter)
# - catch something like f(x,y) for x in domain1, y in domain2
process_generator(gen::Base.Generator{<:Base.Iterators.ProductIterator}) =
    process_generator(gen, gen.iter.iterators)
function process_generator(gen, iterators::Tuple{Vararg{AnyDomain}})
    domain = productdomain(map(domain, iterators)...)
    dims = map(dimension, iterators)
    (gen.f, domain)
end

integrate(integrand, args...) =
    integrate1(integrand, process_arguments(args...)...)

integrate1(integrand, domain, measure, properties...) =
    integrate(suggestedstrategy(domain, measure, properties...), integrand, domain, measure, properties...)

integrate(qs::QuadratureStrategy, integrand, args...) =
    integrate_start(qs, integrand, process_arguments(args...)...)


# The process is as follows:
# - The argument list is completed (missing domain or measure are added)
#   -> integrate_property is invoked
# - integrate_property processes properties using dispatch
#   -> integrate_measure is invoked
# - integrate_measure processes measures using dispatch
#   -> integrate_domain is invoked
# - integrate_domain processes domains using dispatch
#   -> apply_quad is invoked
# - If apply_quad is not intercepted, a fallback routine is invoked.


integrate_start(qs, integrand, domain, measure, properties...) =
    integrate_property(qs, FunIntegrand{promote_type(numtype(DomainRef(domain)),codomaintype(measure))}(integrand), domain, measure, properties...)

integrate_start(qs, integrand::Integrand, domain, measure, properties...) =
    integrate_property(qs, integrand, domain, measure, properties...)


function sum_integrals(qs, integrand::Integrand, domains, measure, properties...)
    Itot, Etot = zero_result(integrand)
    for domain in domains
        I,E = integrate_start(qs, integrand, domain, measure, filter_singularities(properties...)...)
        Itot += I
        Etot += I
    end
    Itot, Etot
end

# Step 1: process any singularities in the list of properties
function integrate_property(qs::AdaptiveStrategy, integrand, domain, measure, property, properties...)
    if property_splits_domain(property, properties...)
        domains = property_split_domain(domain, property, properties...)
        sum_integrals(qs, integrand, domains, measure, filter_singularities(property, properties...)...)
    else
        transformed_arguments = process_properties(qs, integrand, domain, measure, property, properties...)
        integrate_measure(qs, transformed_arguments...)
    end
end

integrate_property(qs, integrand, domain, measure, properties...) =
    integrate_measure(qs, integrand, domain, measure, properties...)

# STEP M: dispatch on the measure
function integrate_measure(qs, integrand, domain, measure, properties...)
    transformed_arguments = process_measure(qs, integrand, domain, measure, properties...)
    integrate_domain(qs, transformed_arguments...)
end

# STEP D: dispatch on the domain
function integrate_domain(qs, integrand, domain, measure, properties...)
    if domain_splits(domain)
        domains = domain_split(domain)
        sum_integrals(qs, integrand, domains, measure, properties...)
    else
        transformed_arguments = process_domain(qs, integrand, domain, measure, properties...)
        if domain_splits(transformed_arguments[2])
            integrate_domain(qs, transformed_arguments...)
        else
            integrate_done(qs, transformed_arguments...)
        end
    end
end

integrate_done(qs, integrand::FunIntegrand, domain, measure, properties...) =
    integrate_done(qs, integrand.fun, domain, measure, properties...)
integrate_done(qs, integrand, domain, measure, properties...) =
    select_quad(qs, integrand, domain, measure, properties...)

# Selection: use a suitable quadrature strategy

# Generic adaptive quadrature in 1D invokes QuadGK, everywhere else it uses hcubature
select_quad(qs::QuadAdaptive, integrand, domain::Domain{T}, measure, properties...) where {T<:Number} =
    apply_quad(Q_quadgk(qs), integrand, domain, measure, properties...)
select_quad(qs::QuadAdaptive, integrand, domain::Domain{T}, measure, properties...) where {T} =
    apply_quad(Q_hcubature(qs), integrand, domain, measure, properties...)
select_quad(qs, integrand, domain, measure, properties...) =
    apply_quad(qs, integrand, domain, measure, properties...)

# FINAL STEP: invoke apply_quad, with a fallback if necessary
apply_quad(qs, integrand, domain, measure, properties...) =
    fallback_integrate(qs, integrand, domain, measure, properties...)

fallback_integrate(qs, integrand, domain, measure, properties...) =
    error("Don't know how to integrate on $domain with strategy $(qs).")


# For quadgk, we only know how to compute intervals
function apply_quad(qs::Q_quadgk, integrand, domain::AbstractInterval, measure::LebesgueMeasure, properties...)
    if isempty(domain)
        zero_result(integrand, prectype(DomainRef(domain)))
    else
        quadgk(integrand, extrema(domain)...; atol = qs.atol, rtol = qs.rtol, maxevals = qs.maxevals)
    end
end

function apply_quad(qs::Q_hcubature, integrand, domain, measure, properties...)
    if islebesguemeasure(measure)
        apply_hcubature(qs, integrand, domain)
    else
        error("HCubature is invoked, but not with a Lebesgue measure: $(measure)")
    end
end

function apply_hcubature(qs::Q_hcubature, integrand, domain)
    f = x -> integrand(x) * DomainSets.indomain(x, domain)
    apply_hcubature(qs, f, boundingbox(domain))
end

function apply_hcubature(qs::Q_hcubature, integrand, domain::ProductDomain)
    if all(DomainSets.isinterval, components(domain))
        a = map(float, map(infimum, components(domain)))
        b = map(float, map(supremum, components(domain)))
        hcubature(integrand, a, b; rtol=qs.rtol, atol=qs.atol, maxevals = qs.maxevals)
    else
        apply_hcubature(qs, x -> integrand(x)*DomainSets.indomain(x, domain), boundingbox(domain))
    end
end

# Given a rule on [-1,1] and a different interval, we map the rule
function apply_quad(qs::ChebyshevIntervalRule, integrand, domain::AbstractInterval, measure::LebesgueMeasure, properties...)
    a = infimum(domain)
    b = supremum(domain)
    D = (b-a)/2
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(a + (x[i]+1)*D)*D for i in 1:length(x))
    z, unknown_error(z)
end

function apply_quad(qs::UnitIntervalRule, integrand, domain::AbstractInterval, measure::LebesgueMeasure, properties...)
    a = leftendpoint(domain)
    b = rightendpoint(domain)
    D = b-a
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(a + x[i]*D)*D for i in 1:length(x))
    z, unknown_error(z)
end

function apply_quad(qs::IntervalRule, integrand, interval::AbstractInterval, measure::LebesgueMeasure, properties...)
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
function apply_quad(qs::HalfLineRule, integrand, domain::HalfLine, measure::LebesgueMeasure, properties...)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end

# Make sure that a real line rule is applied to a real line integral
function apply_quad(qs::RealLineRule, integrand, domain::RealLine, measure::LebesgueMeasure, properties...)
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end

# Any other quadrature rule, we just apply it
function apply_quad(qs::GenericDomainRule, integrand, domain, measure, properties...)
    @assert domain(qs) == domain
    x = points(qs)
    w = weights(qs)
    z = sum(w[i]*integrand(x[i]) for i in 1:length(x))
    z, unknown_error(z)
end
