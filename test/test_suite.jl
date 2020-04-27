
using Test
using DomainSets
using HCubature
using StaticArrays

using DomainIntegrals

include("test_measures.jl")
include("test_integrals.jl")
include("test_gauss.jl")

@testset "Measures" begin
    test_measures()
end

@testset "Integrals" begin
    test_integrals()
end

@testset "Quadrature rules" begin
    test_gauss()
end
