
## Process domains

"Is the domain associated with a simpler domain for numerical integration?"
has_integration_domain(d) = hasparameterization(d)
"Return the integration domain associated with `d`."
integration_domain(d) = parameterdomain(d)
"Return the map from the integration domain of `d` to `d`."
mapfrom_integration_domain(d) = mapfrom_parameterdomain(d)

# Avoid the affine map for intervals and cubes
has_integration_domain(d::Interval) = false
integration_domain(d::Interval) = d
has_integration_domain(d::DomainSets.HyperRectangle) = false
integration_domain(d::DomainSets.HyperRectangle) = d

"Does the integration domain naturally split into subdomains?"
domain_splits(domain) = false
domain_splits(domain::UnionDomain) = true
domain_splits(domain::ProductDomain) = any(map(domain_splits, components(domain)))

"Split the integration domain into subdomains."
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

"Convert a union of domains into a vector of domains without overlap."
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

function process_domain(qs, integrand, domain::DomainSets.TupleProductDomain, measure::Lebesgue, properties...)
    error("Don't know how to integrate on a tuple product domain.")
end

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
#     process_domain(qs, (t->jacdet(m,t)) * (integrand ∘ m), superdomain(domain), measure, properties...)
# end

# enable integration over the unit simplex
has_integration_domain(::EuclideanUnitSimplex{2}) = true
integration_domain(d::EuclideanUnitSimplex{2,T}) where T = UnitSquare{T}()
mapfrom_integration_domain(d::EuclideanUnitSimplex{2,T}) where T =
    DuffyTransform{T}()

# map upper and lower right triangles to the unit simplex
has_integration_domain(::LowerRightTriangle) = true
integration_domain(::LowerRightTriangle{T}) where T = EuclideanUnitSimplex{2,T}()
function mapfrom_integration_domain(d::LowerRightTriangle{T}) where T
    scalefactor = d.b-d.a
    # from the unit simplex to the lower right triangle in the unit square
    m1 = AffineMap(SA{T}[-1 0; 0 1], SA{T}[1; 0])
    # from the unit square to [a,b]^2
    m2 = AffineMap(SA{T}[scalefactor 0; 0 scalefactor], SA{T}[d.a; d.a])
    m2 ∘ m1
end

has_integration_domain(::UpperRightTriangle) = true
integration_domain(::UpperRightTriangle{T}) where T = EuclideanUnitSimplex{2,T}()
function mapfrom_integration_domain(d::UpperRightTriangle{T}) where T
    scalefactor = d.b-d.a
    m1 = AffineMap(SA{T}[1 0; 0 -1], SA{T}[0; 1])
    m2 = AffineMap(SA{T}[scalefactor 0; 0 scalefactor], SA{T}[d.a; d.a])
    m2 ∘ m1
end
