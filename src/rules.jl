
using DomainSets: element, elements, numelements, ×

recombine_outcome(IEs) = sum(IE[1] for IE in IEs), sum(IE[2] for IE in IEs)

## Singularities

splitdomain_sing(sing, domain) = (domain,)

splitdomain_sing(sing::PointSingularity, domain) = splitdomain_point(point(sing), domain)

splitdomain_point(point, domain::Domain) = (domain,)

function splitdomain_point(x::T, domain::AbstractInterval{T}) where {T}
    tol = tolerance(T)
    a, b = extrema(domain)
    if (x > a+tol) && (x < b-tol)
        # x is in [a,b] and sufficiently far away from the boundary
        (a..x-tol, x+tol..b)
    elseif (x < a-tol) || (x > b+tol)
        # x is not in [a,b] and sufficiently far away from the boundary
        (domain,)
    else
        # x is close to a or b
        if abs(x-a) <= tol
            (max(a,x)..b,)
        else
            (a..min(b,x),)
        end
    end
end

splitdomain_point(x, domain::ProductDomain) = splitdomain_point(x, domain, elements(domain)...)

function splitdomain_point(x, domain::ProductDomain, domain1, domain2)
    domains1 = splitdomain_point(x[1], domain1)
    domains2 = splitdomain_point(x[2], domain2)
    [d1 × d2 for d1 in domains1,d2 in domains2]
end

splitdomain_sing(sing::DiagonalSingularity, domain) =
    splitdomain_diagonal(domain)

splitdomain_diagonal(domain::ProductDomain) = splitdomain_diagonal(domain, elements(domain)...)

function splitdomain_diagonal(domain::ProductDomain, domain1::AbstractInterval, domain2::AbstractInterval)
    diff1 = domain1 \ domain2
    diff2 = domain2 \ domain1
    overlap = domain1 ∩ domain2
    if width(overlap) < 100eps(eltype(domain1))
        # "if isempty(overlap)" is not sufficient, because the overlap may be numerically small
        [domain1 × domain2]
    else
        a, b = extrema(overlap)
        if isempty(diff1) && isempty(diff2)
            [LowerRightTriangle(a, b), UpperRightTriangle(a, b)]
        elseif isempty(diff1)
            [LowerRightTriangle(a, b), UpperRightTriangle(a, b), domain1 × diff2]
        elseif isempty(diff2)
            [diff1 × domain2, LowerRightTriangle(a, b), UpperRightTriangle(a, b)]
        else
            [diff1 × domain2, LowerRightTriangle(a, b), UpperRightTriangle(a, b), overlap × diff2]
        end
    end
end


const SplittingSingularity = Union{PointSingularity,DiagonalSingularity}

function quadrature_s(qs, integrand, domain, measure, sing::SplittingSingularity)
    domains = splitdomain_sing(sing, domain)
    if length(domains) > 1
        IEs = quadrature_m.(Ref(qs), Ref(integrand), domains, Ref(measure), Ref(sing))
        recombine_outcome(IEs)
    else
        quadrature_m(qs, integrand, domain, measure, sing)
    end
end




## Measures

function quadrature_m(qs, integrand, domain, δ::DiracMeasure, sing)
    x = point(δ)
    if x ∈ domain
        I = integrand(x)
        E = zero(I)
    else
        I, E = zero_result(integrand, eltype(domain))
    end
    I, E
end

function quadrature_m(qs, integrand, domain, μ::DiscreteMeasure, sing)
    x = points(μ)
    w = weights(μ)
    I,E = zero_result(integrand, eltype(domain))
    for i in 1:length(x)
        if x[i] ∈ domain
            I += w[i]*integrand(x[i])
        end
    end
    I, E
end

# We replace all measures by the LebesgueMeasure on the space
function quadrature_m(qs, integrand, domain, measure::Measure{T}, sing) where {T}
    integrand2 = t -> integrand(t) * unsafe_weight(measure, t)
    quadrature_d(qs, integrand2, domain, LebesgueMeasure{T}(), sing)
end

# Truncate an infinite domain to a finite one for numerical evaluation of Hermite integrals
function quadrature_m(qs::AdaptiveStrategy, integrand, domain::FullSpace{T}, measure::HermiteMeasure{T}, sing) where {T}
    integrand2 = t -> integrand(t) * unsafe_weight(measure, t)
    U = sqrt(-log(eps(Float64)))
    quadrature_d(qs, integrand2, -U..U, LebesgueMeasure{T}(), sing)
end

# apply the cosine map for integrals with the ChebyshevT weight, to avoid the singularities
function quadrature_m(qs::AdaptiveStrategy, integrand, domain::ChebyshevInterval, measure::ChebyshevTMeasure{T}, sing) where {T}
    Tpi = convert(T, pi)
    integrand2 = t -> Tpi*integrand(cos(Tpi*t))
    quadrature_d(qs, integrand2, UnitInterval{T}(), UnitLebesgueMeasure{T}(), sing)
end

# apply the cosine map for integrals with the ChebyshevU weight as well
function quadrature_m(qs::AdaptiveStrategy, integrand, domain::ChebyshevInterval, measure::ChebyshevUMeasure{T}, sing) where {T}
    Tpi = convert(T, pi)
    integrand2 = t -> Tpi*integrand(cos(Tpi*t))*sin(Tpi*t)^2
    quadrature_d(qs, integrand2, UnitInterval{T}(), UnitLebesgueMeasure{T}(), sing)
end

# Lebesgue measures can pass through unaltered
quadrature_m(qs, integrand, domain, measure::AbstractLebesgueMeasure, sing) =
    quadrature_d(qs, integrand, domain, measure, sing)

# The "best" rule for certain measures becomes a Gauss rule
quadrature_m(qs::BestRule, integrand, domain::ChebyshevInterval, μ::LegendreMeasure, sing) =
    quadrature_m(Q_GaussLegendre(qs.n), integrand, domain, μ, sing)

quadrature_m(qs::BestRule, integrand, domain::ChebyshevInterval, μ::JacobiMeasure, sing) =
    quadrature_m(Q_GaussJacobi(qs.n, μ.α, μ.β), integrand, domain, LegendreMeasure(), sing)

quadrature_m(qs::BestRule, integrand, domain::HalfLine, μ::LaguerreMeasure{T}, sing) where {T} =
    quadrature_m(Q_GaussLaguerre(qs.n, μ.α), integrand, domain, LebesgueMeasure{T}(), sing)

quadrature_m(qs::BestRule, integrand, domain::FullSpace{T}, μ::HermiteMeasure{T}, sing) where {T} =
    quadrature_m(Q_GaussHermite(qs.n), integrand, domain, LebesgueMeasure{T}(), sing)

# function quadrature_m(qs, integrand, domain::AbstractInterval, measure::ChebyshevMeasure, sing)
#     a, b = extrema(domain)
#     quadrature(qs, t->integrand(cos(t))*sin(t), acos(a)..acos(b), LegendreMeasure(), sing)
# end


## Domains

quadrature_d(qs, integrand, domain::EmptySpace, measure, sing) =
    zero_result(integrand, eltype(domain))

function nonoverlapping_domains(domain::UnionDomain)
    domains = [element(domain, 1)]
    for i in 2:numelements(domain)
        d = element(domain, i)
        for j=1:i-1
            d = d \ domains[j]
        end
        push!(domains, d)
    end
    domains
end


expand_domain(domain::Domain) = (domain,)

expand_domain(domain::UnionDomain) = nonoverlapping_domains(domain)

expand_domain(domain::ProductDomain) = expand_domain(domain, elements(domain)...)

function expand_domain(domain::ProductDomain, domain1::Domain, domain2::Domain)
    domains1 = expand_domain(domain1)
    domains2 = expand_domain(domain2)
    [d1 × d2 for d1 in domains1,d2 in domains2]
end

function expand_domain(domain::ProductDomain, domain1::Domain, domain2::Domain, domain3::Domain)
    domains1 = expand_domain(domain1)
    domains2 = expand_domain(domain2)
    domains3 = expand_domain(domain3)
    [×(d1, d2, d3) for d1 in domains1,d2 in domains2,d3 in domains3]
end

function expand_domain(domain::ProductDomain, domain1::Domain, domain2::Domain, domain3::Domain, domain4::Domain)
    domains1 = expand_domain(domain1)
    domains2 = expand_domain(domain2)
    domains3 = expand_domain(domain3)
    domains3 = expand_domain(domain3)
    [×(d1,d2,d3,d4) for d1 in domains1,d2 in domains2,d3 in domains3]
end


const ExpandableDomain = Union{UnionDomain,ProductDomain}

function quadrature_d(qs, integrand, domain::ExpandableDomain, measure, sing)
    domains = expand_domain(domain)
    if length(domains) > 1
        # We deliberately invoke quadrature here, starting the chain again from the start.
        # This gives the earlier steps a change to simplify again, with the new domains
        IEs = quadrature.(Ref(qs), Ref(integrand), domains, Ref(measure), Ref(sing))
        recombine_outcome(IEs)
    else
        # Nothing changed, we move on in the chain with the original domain
        select_quad(qs, integrand, domain, measure, sing)
    end
end

function quadrature_d(qs, integrand, domain::LowerRightTriangle{T}, measure::LebesgueMeasure, sing) where {T}
    # Restriction to Lebesgue measure because we would have to map the measure too
    # TODO with a more general framework for change-of-variables in integrals
    # For now, we hardcode the change of variables and we use Duffy's trick.
    a = domain.a
    b = domain.b
    # For x from a to b, and y from a to x, we set y = a+(x-a)*u, where u goes from 0 to 1. We then have dy = (x-1)du.
    d = UnitInterval{T}()
    square = (a..b) × d
    # TODO: change the singularity to something sensible: the diagonal is mapped to the right side of the square
    quadrature(qs, x -> integrand(SVector(x[1],a+x[2]*(x[1]-a))) * (x[1]-a),
        square, measure, NoSingularity())
end

function quadrature_d(qs, integrand, domain::UpperRightTriangle{T}, measure::LebesgueMeasure, sing) where {T}
    a = domain.a
    b = domain.b
    # For x from a to b, and y from x to b, we set y = b-(b-x)*u, where u goes from 0 to 1. We then have dy = (b-x)du.
    d = UnitInterval{T}()
    square = (a..b) × d
    quadrature(qs, x -> integrand(SVector(x[1], b-(b-x[1])*(1-x[2])))*(b-x[1]),
        square, measure, NoSingularity())
end
