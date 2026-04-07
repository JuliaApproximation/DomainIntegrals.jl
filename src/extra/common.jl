
"A triangle defined by three points in a plane."
abstract type Triangle{T} <: Domain{SVector{2,T}} end
# This is just a conceptual type that acts as a data container. No functionality is implemented
# to manipulate triangles.

"The triangle below the diagonal of a square `[a,b] × [a,b]`."
struct LowerRightTriangle{T} <: Triangle{T}
   a  :: T
   b  :: T
end

DomainSets.indomain(x, d::LowerRightTriangle) =
   (x[1] >= d.a) && (x[1] <= d.b) && (x[2] >= d.a) && (x[2] <= d.b) && (x[2] <= x[1])

"The triangle above the diagonal of a square `[a,b] × [a,b]`."
struct UpperRightTriangle{T} <: Triangle{T}
   a  :: T
   b  :: T
end

DomainSets.indomain(x, d::UpperRightTriangle) =
   (x[1] >= d.a) && (x[1] <= d.b) && (x[2] >= d.a) && (x[2] <= d.b) && (x[2] >= x[1])

const StaticTypes = DomainSets.StaticTypes

"""
The Duffy transform maps the unit square to the unit simplex.

It does so by collapsing the rightmost edge onto the x-axis.
"""
struct DuffyTransform{T} <: DomainSets.Map{SVector{2,T}}
end

DuffyTransform() = DuffyTransform{Float64}()

FunctionMaps.similarmap(::DuffyTransform, ::Type{SVector{2,T}}) where T =
   DuffyTransform{T}()

FunctionMaps.applymap(::DuffyTransform{T}, x) where T =
   SVector{2,T}(x[1],(1-x[1])*x[2])

FunctionMaps.jacobian(::DuffyTransform{T}, x) where T = SA{T}[1 0; -x[2] 1-x[1]]
FunctionMaps.jacdet(::DuffyTransform{T}, x) where T = one(T)-x[1]
FunctionMaps.mapsize(::DuffyTransform) = (2,2)
FunctionMaps.diffvolume(::DuffyTransform) = x -> 1-x[1]
