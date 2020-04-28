
"A triangle defined by three points in a plane."
abstract type Triangle end
# This is just a conceptual type that acts as a data container. No functionality is implemented
# to manipulate triangles.

"The triangle below the diagonal of a square `[a,b] × [a,b]`."
struct LowerRightTriangle{T} <: Triangle
   a  :: T
   b  :: T
end

"The triangle above the diagonal of a square `[a,b] × [a,b]`."
struct UpperRightTriangle{T} <: Triangle
   a  :: T
   b  :: T
end
