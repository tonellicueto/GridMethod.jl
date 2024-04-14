module Polynomial
using LinearAlgebra

export PolynomialSystem
export normalizePoly

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

function normalizePoly(P::PolynomialSystem{T}, polyNorm) where T
    invnorms = [
        polyNorm(P.degrees[i], P.monomialDegrees[i], P.coefficients[i])
        for i in eachindex(P.degrees)
    ]
    invnormsDiag = Diagonal(invnorms)
    return PolynomialSystem{T}(
        v -> invnorms.*P(v),
        v -> invnormsDiag*P.jacobian(v),
        P.degrees,
        P.coefficients.*invnorms,
        P.monomialDegrees
    )
end
end