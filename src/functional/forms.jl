
export Integral

"""
A one-form is a map `V → K`, i.e., it is a generic function of one argument.

A one-form is applied to an argument using the `F(f)` syntax, which itself
invokes `form(F::OneForm, f)`. The `form` function can be specialized for
concrete types of `F` and `f`.

The main use cases in `DomainIntegrals.jl` are a lazy `Integral` object and a
quadratic form that derives from a bilinear form.
"""
abstract type OneForm end

(F::OneForm)(f) = form(F, f)


"`Integral` is a functional that represents integration of a function."
abstract type AbstractIntegral <: OneForm end

domain(F::AbstractIntegral) = F.domain
measure(F::AbstractIntegral) = F.measure
properties(F::AbstractIntegral) = F.properties

form(F::AbstractIntegral, f) = integral(f, domain(F), measure(F), properties(F)...)

struct Integral <: AbstractIntegral
    integrand
    domain
    measure
    properties
end


Integral(f, args...) = Integral(f, process_arguments(args...)...)

Integral(integrand, domain::Domain, measure::Measure, properties::Property...) =
    Integral(integrand, domain, measure, properties)


"""
A two-form is a map `V x W → K`, i.e., it is a generic function of two arguments.

The two-form is applied to two arguments using the `F(f,g)` syntax, which itself
invokes `form(F::TwoForm, f, g)`. The `form` function can be specialized for
concrete types of `F`, `f` and `g`.
"""
abstract type TwoForm end

(F::TwoForm)(f, g) = form(F, f, g)


"A bilinear form is linear in its two arguments."
abstract type BilinearForm <: TwoForm end

"A sesquilinear form is linear in the first and conjugate linear in the second argument."
abstract type SesquilinearForm <: TwoForm end

isbilinear(F::BilinearForm) = true
isbilinear(F::SesquilinearForm) = false

issesquilinear(F::BilinearForm) = false
issesquilinear(F::SesquilinearForm) = true


"""
Supertype of quadratic forms that derive from a two-form.

If `F` is a TwoForm, then we call `F(f,f)` the associated quadratic form.
"""
abstract type AbstractQuadraticForm <: OneForm end

# for convenience, assume the underlying form is stored in the field F
superform(F::AbstractQuadraticForm) = F.F

form(F::AbstractQuadraticForm, f) = quadraticform(superform(F), f)
quadraticform(F, f) = form(F, f, f)

"A quadratic form derived from a stored two-form."
struct QuadraticForm{FORM} <: AbstractQuadraticForm
    F   ::  FORM
end



"""
Apply a two-form to the given arguments.

`form(F::TwoForm, f, g)`

By default, `form` invokes `form1` which in turn invokes `form2`. The `form`
function itself can be specialized on the type of `F` without ambiguity.
Similarly, `form1` and `form2` can be specialized without ambiguity on the types
of `f` and `g` respectively.

Each of these functions can be specialized on all three types involved, as long
as the main type is at least as specific as any other specialization.
"""
form(F::TwoForm, f, g) = form1(F, f, g)

"Helper function for `form(F, f, g)` that can be specialized on the type of `f`."
form1(F, f, g) = form2(F, f, g)

"Helper function for `form(F, f, g)` that can be specialized on the type of `g`."
form2(F, f, g) = default_form(F, f, g)

functionproduct(f, g) = x -> f(x)*g(x)
functionconjugateproduct(f, g) = x -> f(x)*conj(g(x))

# By default we assume that the form is induced by a measure.
default_form(F::BilinearForm, f, g) =
    integral(functionproduct(f,g), measure(F))
default_form(F::SesquilinearForm, f, g) =
    integral(functionconjugateproduct(f,g), measure(F))


##################
## Inner products
##################

"An inner-product is a positive definite sesquilinear form."
abstract type InnerProduct <: SesquilinearForm end

# For convenience, assume the underlying form is stored in the field F.
superform(F::InnerProduct) = F.F

ishermitian(F::InnerProduct) = true
ispositivesemidefinite(F::InnerProduct) = true
ispositivedefinite(F::InnerProduct) = true

# an inner product is typically implemented in terms of another form
form(F::InnerProduct, f, g) = form(superform(F), f, g)

# An example concrete type follows. Other concrete types may be defined that
# return a form on-the-fly by defining `superform`, for example.

"A generic inner product that stores a positive definite Hermitian sesquilinear form."
struct GenericInnerProduct{FORM} <: InnerProduct
    F   ::  FORM
end


#########################
## Measure induced forms
#########################

"""
A bilinear form induced by a measure.

It is the case that `F(f,g)` equals `integral(x -> f(x)*g(x), measure(F))`.
"""
abstract type InducedBilinearForm <: BilinearForm end

issymmetric(F::InducedBilinearForm) = true
isconjugatesymmetric(F::InducedBilinearForm) = false
ishermitian(F::InducedBilinearForm) = false


"""
A sesquilinear form induced by a measure.

It is the case that `F(f,g)` equals `integral(x -> f(x)*conj(g(x)), measure(F))`.
"""
abstract type InducedSesquilinearForm <: SesquilinearForm end

issymmetric(F::InducedSesquilinearForm) = false
isconjugatesymmetric(F::InducedSesquilinearForm) = true
ishermitian(F::InducedSesquilinearForm) = isreal(measure(F))
ispositivesemidefinite(F::InducedSesquilinearForm) = ishermitian(F)


# Two concrete types follow. Other concrete types may be defined that do not
# store a measure, but merely return it on-the-fly by implementing `measure(F)`.

"A bilinear form induced by a generic measure."
struct MeasureBilinearForm{M} <: InducedBilinearForm
    measure     ::  M
end

"A sesquilinear form induced by a generic measure."
struct MeasureSesquilinearForm{M} <: InducedSesquilinearForm
    measure     ::  M
end
