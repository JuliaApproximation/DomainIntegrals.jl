
## Process domains

"Associate the domain with a simpler domain for numerical integration."
has_integration_domain(d) = hasparameterization(d)
"Return the integration domain associated `d`."
integration_domain(d) = parameterdomain(d)
"Return the map from the integration domain of `d` to `d`."
mapfrom_integration_domain(d) = mapfrom_parameterdomain(d)

# Avoid the affine map for intervals and cubes
has_integration_domain(d::Interval) = false
has_integration_domain(d::DomainSets.HyperRectangle) = false

domain_splits(domain) = false
domain_splits(domain::UnionDomain) = true
domain_splits(domain::ProductDomain) = any(map(domain_splits, components(domain)))

domain_split(domain) = (domain,)

domain_split(domain::UnionDomain) = nonoverlapping_domains(domain)
domain_split(domain::ProductDomain) = domain_split(domain, components(domain)...)

function domain_split(domain::ProductDomain, domain1::Domain, domain2::Domain)
    domains1 = domain_split(domain1)
    domains2 = domain_split(domain2)
    [productdomain(d1, d2) for d1 in domains1,d2 in domains2]
end

function domain_split(domain::ProductDomain, domain1::Domain, domain2::Domain, domain3::Domain)
    domains1 = domain_split(domain1)
    domains2 = domain_split(domain2)
    domains3 = domain_split(domain3)
    [productdomain(d1, d2, d3) for d1 in domains1, d2 in domains2, d3 in domains3]
end

function domain_split(domain::ProductDomain, domain1::Domain, domain2::Domain, domain3::Domain, domain4::Domain)
    domains1 = domain_split(domain1)
    domains2 = domain_split(domain2)
    domains3 = domain_split(domain3)
    domains4 = domain_split(domain4)
    [productdomain(d1,d2,d3,d4) for d1 in domains1, d2 in domains2, d3 in domains3, d4 in domains4]
end

"Convert a union of domains into a vector of domains without overlap"
function nonoverlapping_domains(domain::UnionDomain{T}) where {T}
    domains = Vector{Domain{T}}()
    push!(domains, component(domain, 1))
    for i in 2:ncomponents(domain)
        d = component(domain, i)
        for j=1:i-1
            d = d \ domains[j]
        end
        push!(domains, d)
    end
    domains
end

process_domain(qs, integrand, domain, measure, properties...) =
    (integrand, domain, measure, properties...)

function process_domain(qs, integrand, domain::Domain, measure::Lebesgue, properties...)
    if has_integration_domain(domain)
        paramdomain = integration_domain(domain)
        fmap = mapfrom_integration_domain(domain)
        process_domain(qs, diffvolume(fmap) * (integrand ∘ fmap), paramdomain, Lebesgue{eltype(paramdomain)}(), properties...)
        # (diffvolume(fmap) * (integrand ∘ fmap), paramdomain, Lebesgue{eltype(paramdomain)}(), properties...)
    else
        (integrand, domain, measure, properties...)
    end
end


integrate_domain(qs, integrand, domain::EmptySpace, measure, properties...) =
    zero_result(integrand)

integrate_domain(qs, integrand, domain::Point, measure, properties...) =
    zero_result(integrand)

# function process_domain(qs, integrand, domain::MappedDomain, measure::Lebesgue, properties...)
#     m = forward_map(domain)
#     @show m
#     process_domain(qs, (t->jacdet(m,t)) * (integrand ∘ m), superdomain(domain), measure, properties...)
# end

# TODO: implement triangles in terms of UnitSimplex, rather than the other way around
function process_domain(qs, integrand, domain::EuclideanUnitSimplex{2,T}, measure::Lebesgue, properties...) where {T}
    m = AffineMap(SA[-1 zero(T); 0 1], SA[one(T); 0])
    process_domain(qs, diffvolume(m) * (integrand ∘ m), LowerRightTriangle(zero(T),one(T)), measure, properties...)
end

function process_domain(qs, integrand, domain::LowerRightTriangle{T}, measure::Lebesgue, properties...) where {T}
    a = domain.a
    b = domain.b
    # For x from a to b, and y from a to x, we set y = a+(x-a)*u, where u goes from 0 to 1. We then have dy = (x-1)du.
    d = UnitInterval{T}()
    square = (a..b) × d
    triangle_map = x -> SVector(x[1],a+x[2]*(x[1]-a))
    # TODO: change the singularity to something sensible: the diagonal is mapped to the right side of the square
    ((x->(x[1]-a)) * (integrand ∘ triangle_map), square, measure)
end

function process_domain(qs, integrand, domain::UpperRightTriangle{T}, measure::Lebesgue, properties...) where {T}
    a = domain.a
    b = domain.b
    # For x from a to b, and y from x to b, we set y = b-(b-x)*u, where u goes from 0 to 1. We then have dy = (b-x)du.
    d = UnitInterval{T}()
    square = (a..b) × d
    triangle_map = x -> SVector(x[1], b-(b-x[1])*(1-x[2]))
    ((x->(b-x[1])) * (integrand ∘ triangle_map), square, measure)
end
