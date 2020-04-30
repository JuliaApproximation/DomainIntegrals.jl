
export SingularPoint,
    LogPointSingularity,
    AlgebraicPointSingularity,
    DiagonalSingularity

"The supertype of all kinds of singularities of an integrand."
abstract type Singularity end

"The integrand has no (known) singularity."
struct NoSingularity <: Singularity end

"Supertype of all kinds of point singularities."
abstract type PointSingularity <: Singularity end

point(s::PointSingularity) = s.point

"An unspecified singularity occurring at a specific point."
struct SingularPoint{T} <: PointSingularity
    point   ::  T
end

"A logarithmic singularity at a specific point."
struct LogPointSingularity{T} <: PointSingularity
    point   ::  T
end

"An algebraic singularity at a specific point."
struct AlgebraicPointSingularity{O,T} <: PointSingularity
    point   ::  T
    order   ::  O
end

order(s::AlgebraicPointSingularity) = s.order

"Supertype of singularities along a curve."
abstract type CurveSingularity <: Singularity end

"For 2D integrands, a singularity along the line `x=y`."
struct DiagonalSingularity <: CurveSingularity
end


"""
A list of breakpoints of the integrand, for example because the integrand
is piecewise smooth.
"""
struct Waypoints{P} <: Singularity
    points  ::  P
end
