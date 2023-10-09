module DomainIntegrals

using QuadGK, LinearAlgebra
using GaussQuadrature, FastGaussQuadrature, HCubature

using StaticArrays
using IntervalSets, DomainSets
using CompositeTypes.Display

import Base: convert, checkbounds, length,
    *, ∘, ==

import DomainSets:
    domaintype,
    codomaintype,
    prectype,
    numtype,
    forward_map,
    components,
    factors

using DomainSets: ×,
    EmptySpace, FullSpace,
    domain_prectype, domain_numtype,
    euclideandimension


include("extra/common.jl")

include("gauss.jl")
include("functional/measure.jl")
include("weights.jl")
include("singularity.jl")
include("strategy.jl")
include("integrand.jl")
include("processing/properties.jl")
include("processing/measure.jl")
include("processing/domain.jl")
include("integral.jl")
include("complexplane.jl")

include("functional/forms.jl")

end # module
