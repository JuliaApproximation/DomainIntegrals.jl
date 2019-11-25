# DomainIntegrals.jl

[![Build Status](https://travis-ci.org/JuliaApproximation/DomainIntegrals.jl.svg?branch=master)](https://travis-ci.org/JuliaApproximation/DomainIntegral.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/gc75y13g0kerxll8?svg=true)](https://ci.appveyor.com/project/dlfivefifty/domainintegrals-jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaApproximation/DomainIntegrals.jl/badge.svg)](https://coveralls.io/github/JuliaApproximation/DomainIntegrals.jl)


DomainIntegrals.jl is a package designed to numerically evaluate integrals on
domains like they are defined by DomainSets.jl.

The package does not include new methods for numerical integration. It relies
on other packages such as QuadGK and HCubature. The methods of those packages
are leveraged to evaluate integrals on more general domains than intervals and
boxes.


## Examples

```julia
julia> using DomainSets, DomainIntegrals

julia> integral(cos, 0..1.0)
0.8414709848078965
```
