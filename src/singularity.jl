
export SingPoint,
    LogSingPoint,
    AlgebraicSingPoint,
    SingularDiagonal

@deprecate DiagonallySingular SingularDiagonal

"""
A `Property` type describes a property of the integrand that is independent
from the other aspects of the integral (domain, measure).

The main examples are singularities of the integrand function, i.e., smoothness
properties.
"""
abstract type Property end

"No property of the integrand is known."
struct NoProperty <: Property end


"The supertype of all kinds of singularities of an integrand."
abstract type Singularity <: Property end

"No singularity of the integrand is known."
struct NoSingularity <: Singularity end

"Supertype of all kinds of point singularities."
abstract type PointSingularity <: Singularity end

point(s::PointSingularity) = s.point

"An unspecified singularity occurring at a specific point."
struct SingPoint{T} <: PointSingularity
    point   ::  T
end

"An unspecified singularity occurring at a corner point."
struct SingCorner{T} <: PointSingularity
    point   ::  T
end

"A logarithmic singularity at a specific point."
struct LogSingPoint{T} <: PointSingularity
    point   ::  T
end

"An algebraic singularity at a specific point."
struct AlgebraicSingPoint{O,T} <: PointSingularity
    point   ::  T
    order   ::  O
end

order(s::AlgebraicSingPoint) = s.order

"Supertype of singularities along a curve."
abstract type CurveSingularity <: Singularity end

"For 2D integrands, a singularity along the line `x=y`."
struct SingularDiagonal <: CurveSingularity
end


"""
A list of breakpoints of the integrand, for example because the integrand
is piecewise smooth.
"""
struct Waypoints{P} <: Property
    points  ::  P
end

points(P::Waypoints) = P.points
