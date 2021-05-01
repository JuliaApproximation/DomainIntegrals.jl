
"""
Transform the integral problem based on the supplied properties.

For example, if the integrand is known to be singular, then the domain may be
split into several components such that there are no singularities in the
interior.
"""
process_properties(qs, integrand, domain, measure, properties...) =
    (qs, integrand, domain, measure, properties...)

#################
## Singularities
#################

filter_singularities() = ()
filter_singularities(property::Singularity) = ()
filter_singularities(property::Property) = (property,)
filter_singularities(property::Singularity, properties...) =
    filter_singularities(properties...)
filter_singularities(property::Property, properties...) =
    (property, filter_singularities(properties...)...)

property_splits_domain(property) = false
property_splits_domain(property::PointSingularity) = true
property_splits_domain(property::SingularDiagonal) = true

property_splits_domain() = false
property_splits_domain(property1, property2, properties...) =
    property_splits_domain(property1) || property_splits_domain(property2, properties...)

property_split_domain(domain) = (domain,)
property_split_domain(domain, property) = splitdomain_sing(property, domain)

function property_split_domain(domain, property, properties...)
    domains = Domain[]
    push!(domains, domain)
    property_split_multiple(domains, property, properties...)
end

property_split_multiple(domains) = domains
function property_split_multiple(domains_to_split, property, properties...)
    if property_splits_domain(property)
        domains = Domain[]
        for domain_to_split in domains_to_split
            for d in property_split_domain(domain_to_split, property)
                push!(domains, d)
            end
        end
        property_split_multiple(domains, properties...)
    else
        property_split_multiple(domains_to_split, properties...)
    end
end



"""
Split the domain into a list of domains, such that any singularity lies
on the boundary of one or more of the constituting parts.

Care is taken to avoid a singularity in the interior due to roundoff errors.
"""
splitdomain_sing(sing, domain) = (domain,)

splitdomain_sing(sing::PointSingularity, domain) = splitdomain_point(point(sing), domain)

function splitdomain_point(point, domain::MappedDomain)
    m = inverse_map(domain)
    map_domain.(Ref(m), splitdomain_point(m(point), superdomain(domain)))
end

function splitdomain_point(point, domain::Domain)
    if point ∈ domain
        @warn("Don't know how to split domain $(domain) at point $(point)")
        (domain,)
    else
        (domain,)
    end
end

"Split the domain according to a singularity at the point `x`."
function splitdomain_point(x, domain::AbstractInterval)
    T = promote_type(typeof(x),prectype(domain))
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
            (max(a,x)+tol..b,)
        else
            (a..min(b,x)-tol,)
        end
    end
end

splitdomain_point(x, domain::ProductDomain) = splitdomain_point(x, domain, components(domain)...)

function splitdomain_point(x, domain::ProductDomain, domain1, domain2)
    domains1 = splitdomain_point(x[1], domain1)
    domains2 = splitdomain_point(x[2], domain2)
    [d1 × d2 for d1 in domains1,d2 in domains2]
end

splitdomain_sing(sing::SingularDiagonal, domain) =
    splitdomain_diagonal(domain)

splitdomain_diagonal(domain::ProductDomain) = splitdomain_diagonal(domain, components(domain)...)

function splitdomain_diagonal(domain::ProductDomain, domain1::AbstractInterval, domain2::AbstractInterval)
    T = prectype(domain)
    diff1 = domain1 \ domain2
    diff2 = domain2 \ domain1
    overlap = domain1 ∩ domain2
    if width(overlap) < 100eps(T)
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
