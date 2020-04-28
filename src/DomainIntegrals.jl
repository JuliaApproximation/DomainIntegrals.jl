module DomainIntegrals

using QuadGK, LinearAlgebra
using GaussQuadrature, FastGaussQuadrature, HCubature

using StaticArrays
using IntervalSets, DomainSets

import Base: convert

import DomainSets:
    EmptySpace,
    FullSpace,
    domaintype,
    codomaintype,
    prectype,
    numtype


include("common.jl")
include("gauss.jl")
include("measure.jl")
include("singularity.jl")
include("strategy.jl")
include("rules.jl")
include("integral.jl")

end # module
