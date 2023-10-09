
## Process domains

using DomainSets:
	simplifies,
	simplify

"A canonical domain for the purposes of evaluating integrals."
struct CanonicalInt <: DomainSets.CanonicalType
end

import DomainSets:
    canonicaldomain,
    mapfrom_canonical,
	mapto_canonical

"Is the domain associated with a simpler domain for numerical integration?"
has_integration_domain(d) = hascanonicaldomain(CanonicalInt(), d)

"Return the integration domain associated with `d`."
integration_domain(d) = canonicaldomain(CanonicalInt(), d)

"Return the map from the integration domain of `d` to `d`."
mapfrom_integration_domain(d) = mapfrom_canonical(CanonicalInt(), d)
mapto_integration_domain(d) = mapto_canonical(CanonicalInt(), d)


# by default we use the parameter domain
canonicaldomain(::CanonicalInt, d) = parameterdomain(d)
mapfrom_canonical(::CanonicalInt, d) = mapfrom_parameterdomain(d)
# by default one only has to implement mapfrom
mapto_canonical(::CanonicalInt, d) = leftinverse(mapfrom_canonical(CanonicalInt(), d))

# we intercept product domains, so that CanonicalInt is passed on to its factors
canonicaldomain(::CanonicalInt, d::ProductDomain) =
	any(map(has_integration_domain, factors(d))) ?
		productdomain(map(integration_domain, factors(d))...) :
		d
mapto_canonical(::CanonicalInt, d::ProductDomain) =
	DomainSets.matching_product_map(d, map(mapto_integration_domain, factors(d)))
mapfrom_canonical(::CanonicalInt, d::ProductDomain) =
	DomainSets.matching_product_map(d, map(mapfrom_integration_domain, factors(d)))

# same for mapped domains
canonicaldomain(ctype::CanonicalInt, d::DomainSets.AbstractMappedDomain) =
	canonicaldomain(ctype, superdomain(d))
mapfrom_canonical(ctype::CanonicalInt, d::DomainSets.AbstractMappedDomain) =
	forward_map(d) ∘ mapfrom_canonical(ctype, superdomain(d))
mapto_canonical(ctype::CanonicalInt, d::DomainSets.AbstractMappedDomain) =
	mapto_canonical(ctype, superdomain(d)) ∘ inverse_map(d)


# Avoid the affine map for intervals and cubes
canonicaldomain(::CanonicalInt, d::Interval) = d
canonicaldomain(::CanonicalInt, d::DomainSets.HyperRectangle) = d
mapfrom_canonical(::CanonicalInt, d::Interval) = IdentityMap{eltype(d)}()
mapto_canonical(::CanonicalInt, d::Interval) = IdentityMap{eltype(d)}()
mapfrom_canonical(::CanonicalInt, d::DomainSets.HyperRectangle) = IdentityMap{eltype(d)}()
mapto_canonical(::CanonicalInt, d::DomainSets.HyperRectangle) = IdentityMap{eltype(d)}()


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

function process_domain(qs, integrand, domain, measure::Lebesgue, properties...)
	# first simplify the domain
	if simplifies(domain)
		process_domain(qs, integrand, simplify(domain), measure, properties...)
	else
		# next find a parametric description of the domain
	    if has_integration_domain(domain)
	        paramdomain = integration_domain(domain)
	        fmap = mapfrom_integration_domain(domain)
	        process_domain(qs, diffvolume(fmap) * (integrand ∘ fmap), paramdomain, Lebesgue{domaineltype(paramdomain)}(), properties...)
	    else
	        (integrand, domain, measure, properties...)
	    end
	end
end


integrate_domain(qs, integrand, domain::EmptySpace, measure, properties...) =
    zero_result(integrand)

integrate_domain(qs, integrand, domain::Point, measure, properties...) =
    zero_result(integrand)


# enable integration over the unit simplex
canonicaldomain(::CanonicalInt, d::EuclideanUnitSimplex{2,T}) where T = UnitSquare{T}()
mapfrom_canonical(::CanonicalInt, d::EuclideanUnitSimplex{2,T}) where T =
    DuffyTransform{T}()


# map upper and lower right triangles to the unit simplex
canonicaldomain(::CanonicalInt, ::LowerRightTriangle{T}) where T = EuclideanUnitSimplex{2,T}()
function mapfrom_canonical(::CanonicalInt, d::LowerRightTriangle{T}) where T
    scalefactor = d.b-d.a
    # from the unit simplex to the lower right triangle in the unit square
    m1 = AffineMap(SA{T}[-1 0; 0 1], SA{T}[1; 0])
    # from the unit square to [a,b]^2
    m2 = AffineMap(SA{T}[scalefactor 0; 0 scalefactor], SA{T}[d.a; d.a])
    m2 ∘ m1
end

canonicaldomain(::CanonicalInt, ::UpperRightTriangle{T}) where T = EuclideanUnitSimplex{2,T}()
function mapfrom_canonical(::CanonicalInt, d::UpperRightTriangle{T}) where T
    scalefactor = d.b-d.a
    m1 = AffineMap(SA{T}[1 0; 0 -1], SA{T}[0; 1])
    m2 = AffineMap(SA{T}[scalefactor 0; 0 scalefactor], SA{T}[d.a; d.a])
    m2 ∘ m1
end
