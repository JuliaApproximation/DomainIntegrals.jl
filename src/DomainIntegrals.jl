module DomainIntegrals

using QuadGK, LinearAlgebra
using GaussQuadrature, FastGaussQuadrature, HCubature

using StaticArrays
using IntervalSets, DomainSets

import Base: convert, checkbounds, length

import DomainSets:
    EmptySpace,
    FullSpace,
    domaintype,
    codomaintype,
    prectype,
    numtype,
    domain


include("common.jl")
include("gauss.jl")
include("measure.jl")
include("weights.jl")
include("singularity.jl")
include("strategy.jl")
include("rules.jl")
include("integral.jl")
include("forms.jl")

end # module
