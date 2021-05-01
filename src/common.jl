
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
   (x[1] >= d.a) && (x[1] <= d.b) && (x[2] >= d.a) && (x[2] <= d.b) && (x[1] <= x[2])

"The triangle above the diagonal of a square `[a,b] × [a,b]`."
struct UpperRightTriangle{T} <: Triangle{T}
   a  :: T
   b  :: T
end

DomainSets.indomain(x, d::UpperRightTriangle) =
   (x[1] >= d.a) && (x[1] <= d.b) && (x[2] >= d.a) && (x[2] <= d.b) && (x[1] >= x[2])
