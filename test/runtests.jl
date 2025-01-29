
using Test
using DomainSets, FunctionMaps, HCubature, LinearAlgebra
using StaticArrays

using DomainIntegrals

include("test_measures.jl")
include("test_integrals.jl")
include("test_gauss.jl")

@testset "Measures" begin
    test_measures()
    test_discrete_measures()
end

@testset "Quadrature rules" begin
    test_gauss()
end

@testset "Selected integrals" begin
    test_integrals()
end
