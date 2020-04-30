# DomainIntegrals.jl

[![Build Status](https://travis-ci.org/JuliaApproximation/DomainIntegrals.jl.svg?branch=master)](https://travis-ci.org/JuliaApproximation/DomainIntegral.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/gc75y13g0kerxll8?svg=true)](https://ci.appveyor.com/project/dlfivefifty/domainintegrals-jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaApproximation/DomainIntegrals.jl/badge.svg)](https://coveralls.io/github/JuliaApproximation/DomainIntegrals.jl)


DomainIntegrals is a package designed to numerically evaluate integrals on
domains like they are defined by the [DomainSets](https://github.com/JuliaApproximation/DomainSets.jl) package.

The package does not include new methods for numerical integration. It relies
on other Julia packages such as [QuadGK](https://github.com/JuliaMath/QuadGK.jl) and [HCubature](https://github.com/JuliaMath/HCubature.jl). The methods of those packages
are leveraged to evaluate integrals on more general domains than intervals and
boxes.


## Examples

Evaluate the integral of `cos` on the interval `[0,1]` using `integral` or `quadrature`. The `integral` function simply returns a value, while `quadrature`
returns both the value and an estimated accuracy (as returned by the underlying packages):
```julia
julia> using DomainSets, DomainIntegrals

julia> integral(cos, 0..1.0)
0.8414709848078965

julia> integral(x -> exp(x[1]+x[2]), (0..1.0)^2)
2.9524924420120535

julia> quadrature(cos, UnionDomain(0..1, 2..3))
(0.07329356604208204, 1.1102230246251565e-16)
```

It is possible to specify singularities of the integrand. The integration domain is split such that the singularity lies on the boundary:
```julia
julia> integral(t -> sin(log(abs(t))), -1..1, LogPointSingularity(0.0))
-1.0000000021051316

julia> using DomainSets: ×

julia> julia> integral( x -> exp(log(abs(x[1]-x[2]))), (2..3) × (1..4), DiagonalSingularity())
2.333333333333333
```

Weighted integrals are supported through the definition of measures. A few standard weight functions are included, in particular those associated with the classical orthogonal polynomials (Legendre, Chebyshev, Jacobi, Laguerre and Hermite):
```julia
julia> integral(cos, ChebyshevTMeasure())
2.403939430634413

julia> integral(t -> cos(t)*1/sqrt(1-t^2), -1.0..1.0)
2.403939410869398
```
For the particular example of the ChebyshevT measure (associated with Chebyshev polynomials of the first kind), the typical cosine map is applied which removes the algebraic endpoint singularities of the weight function, before it is evaluated numerically.

Optionally, as a first argument to `integral` or `quadrature` the user can specify a quadrature strategy. The default is `AdaptiveStrategy`. Explicitly providing this argument allows setting optional parameters:
```julia
julia> I, E = quadrature(QuadAdaptive(atol=1e-3, rtol = 1e-3), t->cos(t^2), 0..10)
(0.6011251848111901, 0.0004364150560137517)
```

A few well-known quadrature rules are included, as provided by the [GaussQuadrature](https://github.com/billmclean/GaussQuadrature.jl) and [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) packages. They have corresponding strategies. For example, the application of a 10-point Gauss-Laguerre rule:
```julia
julia> integral(Q_GaussLaguerre(10), cos)
0.5000005097999486

julia> integral(t -> cos(t)*exp(-t), HalfLine())
0.5
```


The DomainIntegrals package is extensible. The quadrature routine invokes a series of functions (`quadrature_s`, `quadrature_m`, `quadrature_d`) that allow to
dispatch on the type of singularity, measure and domain respectively. The user
can add methods to these functions to teach DomainIntegrals how to evaluate new kinds of integrals. As an example of a rule that is included, the `quadrature_d` function dispatches on the `DomainUnion` type and recursively evaluates the integrals on each of the composing parts separately (if they do not overlap). The cosine map of Chebyshev measures is implemented by specializing `quadrature_m` for the case of a `ChebyshevTMeasure`. See the file `rules.jl` for other examples.
