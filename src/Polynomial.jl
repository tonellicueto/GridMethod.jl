module Polynomial
using LinearAlgebra

export PolynomialSystem

struct PolynomialSystem{T <: Number}
    evaluate::Function
    jacobian::Function
    degrees::Vector{Int}
    coefficients::Vector{Vector{T}}
    monomialDegrees::Vector{Vector{Vector{Integer}}}
end

function (p::PolynomialSystem{T})(x::Vector{T}) where T
    p.evaluate(x)
end

function normalize(P::PolynomialSystem{T}; polyNorm=nothing) where T
    if isnothing(polyNorm)
        polyNorm = v -> norm(v, 1)
    end

    invnorms = [polyNorm(v) for v in P.coefficients]
    invnormsDiag = Diagonal(invnorms)
    return PolynomialSystem{T}(
        v -> invnorms.*P(v),
        v -> invnormsDiag*P.jacobian(v),
        P.degrees,
        P.coefficients.*invnorms
    )
end
end