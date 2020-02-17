
"What is the floating point precision type associated with the given argument type."
prectype(x) = prectype(typeof(x))
prectype(::Type{T}) where {T} = T
prectype(::Type{T}) where {T <: AbstractFloat} = T
prectype(::Type{Complex{T}}) where {T} = prectype(T)
prectype(::Type{<:AbstractArray{T}}) where {T} = prectype(T)

# Estimate the return type of a scalar valued function with arguments
# of type T. This is often useful as a default.
codomaintype(x) = codomaintype(typeof(x))
codomaintype(::Type{T}) where {T} = T
codomaintype(::Type{T}) where {T <: Number} = T
codomaintype(::Type{<:AbstractArray{T}}) where {T} = codomaintype(T)

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
