
gausslegendre(n::Int) = gausslegendre(Float64, n)
gausslegendre(::Type{Float64}, n) = FastGaussQuadrature.gausslegendre(n)
gausslegendre(::Type{T}, n) where {T} = GaussQuadrature.legendre(T, n)

gausslaguerre(n::Int) = gausslaguerre(Float64, n)
gausslaguerre(n::Int, α::T) where {T<:AbstractFloat} = gausslaguerre(T, n, α)
gausslaguerre(n::Int, α) = gausslaguerre(n, float(α))
gausslaguerre(::Type{T}, n::Int) where {T} = gausslaguerre(T, n, zero(T))
gausslaguerre(::Type{Float64}, n, α) = FastGaussQuadrature.gausslaguerre(n, α)

gaussjacobi(n::Int) = gaussjacobi(Float64, n)
gaussjacobi(n::Int, α::Number, β::Number) = gaussjacobi(n, promote(α, β)...)
gaussjacobi(n::Int, α::N, β::N) where {N<:Number} = gaussjacobi(n, float(α), float(β))
gaussjacobi(n::Int, α::T, β::T) where {T<:AbstractFloat} =
    gaussjacobi(T, n, α, β)
gaussjacobi(::Type{T}, n::Int) where {T} = gaussjacobi(T, n, zero(T), zero(T))
gaussjacobi(::Type{Float64}, n, α, β) = FastGaussQuadrature.gaussjacobi(n, α, β)

gausshermite(n::Int) = gausshermite(Float64, n)
gausshermite(::Type{Float64}, n) = FastGaussQuadrature.gausshermite(n)
