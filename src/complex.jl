
cauchy_contour(a, radius = 1, ::Type{T} = real(typeof(a))) where {T} =
    radius*ComplexUnitCircle{T}() .+ a

export cauchy_integral
"""
    cauchy_integral(f [; a = 0, n  = 1, radius = 1])

Evaluate the integral of `n!/(2πi) f(z)/(z-a)^n` along a circular contour
with radius `radius` around the point `a`.

If the integrand is analytic inside the contour except at the pole around `z=a`
then its value should correspond to the derivative `f^(n)(a)`.
"""
function cauchy_integral(f; a = zero(ComplexF64), radius = 1, n = 1)
    T = real(typeof(a))
    factorial(n)/(2*T(pi)*im) *
        integral( z-> f(z)/(z-a)^(n+1), cauchy_contour(a, radius))
end
