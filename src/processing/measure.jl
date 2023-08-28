
## Process measures

integrate_measure(qs, integrand, domain, μ::LebesgueDomain{T}, properties...) where T =
    integrate_measure(qs, integrand, domain ∩ support(μ), Lebesgue{T}(), properties...)

function integrate_measure(qs, integrand, domain, δ::DiracWeight, properties...)
    x = point(δ)
    if x ∈ domain
        I = integrand(x)
        E = zero(I)
    else
        I, E = zero_result(integrand)
    end
    I, E
end


"A dummy domain to be associated with a discrete measure, for internal purposes."
struct DummyDiscreteDomain{T} <: Domain{T} end

function integrate_measure(qs, integrand, domain::DummyDiscreteDomain, μ::DiscreteWeight, properties...)
    I = sum(w*integrand(x) for (w,x) in zip(weights(μ), points(μ)))
    I, zero_error(I)
end

function integrate_measure(qs, integrand, domain, μ::DiscreteWeight, properties...)
    if domain == covering(μ)
        I = sum(w*integrand(x) for (w,x) in zip(weights(μ), points(μ)))
        I, zero_error(I)
    else
        # Only sum over the elements that lie in the domain
        I, E = zero_result(integrand)
        x = points(μ)
        w = weights(μ)
        for i in 1:length(x)
            if x[i] ∈ domain
                I += w[i]*integrand(x[i])
            end
        end
        I, zero_error(I)
    end
end

# Lebesgue measures can pass through unaltered
integrate_measure(qs, integrand, domain, measure::LebesgueMeasure, properties...) =
    integrate_domain(qs, integrand, domain, measure, properties...)

process_measure(qs, integrand, domain, measure::LebesgueMeasure, properties...) =
    (integrand, domain, measure, properties...)

process_measure(qs, integrand, domain, measure::Measure, properties...) =
    process_measure_default(qs, integrand, domain, measure, properties...)

# By default we replace all measures by the Lebesgue on the space
process_measure_default(qs, integrand, domain, measure::Measure{T}, properties...) where {T} =
    ((t->unsafe_weightfun(measure, t)) * integrand, domain, Lebesgue{T}(), properties...)

function process_measure(qs, integrand, domain::MappedDomain, μ::MappedWeight, properties...)
    map1 = forward_map(μ)
    map2 = forward_map(domain)
    if map1 == map2
        # f -> f(m1(x))
        process_measure(qs, integrand ∘ map1, superdomain(domain), supermeasure(μ), properties...)
    else
        process_measure_default(qs, integrand, domain, μ, properties...)
    end
end

function process_measure(qs, integrand, domain, μ::MappedWeight, properties...)
    # f -> f(m(x))
    map1 = forward_map(μ)
    process_measure(qs, integrand ∘ map1, mapped_domain(map1, domain), supermeasure(μ), properties...)
end


# Truncate an infinite domain to a finite one for numerical evaluation of Hermite integrals
function process_measure(qs::AdaptiveStrategy, integrand, domain::Union{RealLine{T},FullSpace{T}}, measure::HermiteWeight{T}, properties...) where {T}
    U = sqrt(-log(eps(T)))
    (hermite_weightfun * integrand, -U..U, Lebesgue{T}(), properties...)
end

# apply the cosine map for integrals with the ChebyshevT weight, to avoid the singularities
function process_measure(qs::AdaptiveStrategy, integrand, domain::ChebyshevInterval, measure::ChebyshevTWeight{T}, properties...) where {T}
    # Transformation is: f(t) -> pi*f(cos(pi*t))
    Tpi = convert(T, pi)
    map = t -> cos(Tpi*t)
    (Tpi * (integrand ∘ map), UnitInterval{T}(), LebesgueUnit{T}(), properties...)
end

# same as above, but on a subinterval
function process_measure(qs::AdaptiveStrategy, integrand, domain::AbstractInterval, measure::ChebyshevTWeight{T}, properties...) where {T}
    Tpi = convert(T, pi)
    map = t -> cos(Tpi*t)
    a, b = extrema(domain)
    a < -1.001 && throw(BoundsError(measure, a))
    b > 1.001 && throw(BoundsError(measure, b))
    # Set a and b to be within [-1,1] in order to avoid errors with acos below
    a = max(a, -1)
    b = min(b, 1)
    (Tpi * (integrand ∘ map), acos(b)/pi..acos(a)/pi, Lebesgue{T}(), properties...)
end

# apply the cosine map for integrals with the ChebyshevU weight as well
function process_measure(qs::AdaptiveStrategy, integrand, domain::ChebyshevInterval, measure::ChebyshevUWeight{T}, properties...) where {T}
    Tpi = convert(T, pi)
    prefactor = t -> Tpi * sin(Tpi*t)^2
    map = t -> cos(Tpi*t)
    (prefactor * (integrand ∘ map), UnitInterval{T}(), LebesgueUnit{T}(), properties...)
end

function process_measure(qs::AdaptiveStrategy, integrand, domain::AbstractInterval, measure::ChebyshevUWeight{T}, properties...) where {T}
    Tpi = convert(T, pi)
    prefactor = t -> Tpi * sin(Tpi*t)^2
    map = t -> cos(Tpi*t)
    a, b = extrema(domain)
    a < -1.001 && throw(BoundsError(measure, a))
    b > 1.001 && throw(BoundsError(measure, b))
    # Set a and b to be within [-1,1] in order to avoid errors with acos below
    a = max(a, -1)
    b = min(b, 1)
    (prefactor * (integrand ∘ map), acos(b)/pi..acos(a)/pi, Lebesgue{T}(), properties...)
end

function measure_fetch_transformations(qs, domain, weight, properties...)
    T = promote_type(numtype(AsDomain(domain)),numtype(weight))
    I = TransformationLogIntegrand{T}()
    I_trans, domain_trans, weight_trans = process_measure(qs, I, domain, weight)
    I_trans.prefactor, I_trans.fun, domain_trans, weight_trans
end

# For product domains, process the measures dimension per dimension
function process_measure(qs::AdaptiveStrategy, integrand, domain::ProductDomain, μ::ProductWeight, properties...)
    if ncomponents(domain) == ncomponents(μ)
        if ncomponents(domain) == 2
            pre1, map1, domain1, μ1 = measure_fetch_transformations(qs, component(domain, 1), component(μ, 1), properties...)
            pre2, map2, domain2, μ2 = measure_fetch_transformations(qs, component(domain, 2), component(μ, 2), properties...)
            prefactor = t -> pre1(t[1])*pre2(t[2])
            map = t -> SA[map1(t[1]), map2(t[2])]
            prefactor * (integrand ∘ map), productdomain(domain1, domain2), productmeasure(μ1, μ2)
        elseif ncomponents(domain) == 3
            pre1, map1, domain1, μ1 = measure_fetch_transformations(qs, component(domain, 1), component(μ, 1), properties...)
            pre2, map2, domain2, μ2 = measure_fetch_transformations(qs, component(domain, 2), component(μ, 2), properties...)
            pre3, map3, domain3, μ3 = measure_fetch_transformations(qs, component(domain, 3), component(μ, 3), properties...)
            prefactor = t -> pre1(t[1])*pre2(t[2])*pre3(t[3])
            map = t -> SA[map1(t[1]), map2(t[2]), map3(t[3])]
            prefactor * (integrand ∘ map), productdomain(domain1, domain2, domain3), productmeasure(μ1, μ2, μ3)
        else
            process_measure_default(qs, integrand, domain, μ, properties...)
        end
    else
        process_measure_default(qs, integrand, domain, μ, properties...)
    end
end


# The "best" rule for certain measures becomes a Gauss rule
integrate_measure(qs::BestRule, integrand, domain::ChebyshevInterval{T}, μ::LegendreWeight, properties...) where {T} =
    integrate_measure(Q_GaussLegendre(T, qs.n), integrand, domain, μ, properties...)

integrate_measure(qs::BestRule, integrand, domain::ChebyshevInterval, μ::JacobiWeight{T}, properties...) where {T} =
    integrate_measure(Q_GaussJacobi(qs.n, μ.α, μ.β), integrand, domain, LegendreWeight{T}(), properties...)

integrate_measure(qs::BestRule, integrand, domain::HalfLine, μ::LaguerreWeight{T}, properties...) where {T} =
    integrate_measure(Q_GaussLaguerre(qs.n, μ.α), integrand, domain, Lebesgue{T}(), properties...)

integrate_measure(qs::BestRule, integrand, domain::Union{FullSpace{T},RealLine{T}}, μ::HermiteWeight{T}, properties...) where {T} =
    integrate_measure(Q_GaussHermite(qs.n), integrand, domain, Lebesgue{T}(), properties...)
