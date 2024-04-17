module NormsPolynomials
using .Iterators
using LinearAlgebra
using HomotopyContinuation.ModelKit #More efficient to use the package that we will be using anyway
const HCMK = ModelKit
using Logging

export PolynomialSystem
export normalizePoly
export PolyNorm1
export hinfNorm
export matrixInfPNorm
export PolyNormW
export GridSystem

#Norms

function hinfNorm(x::Vector{T}) where T <: Number
    return maximum((one(T), norm(x, Inf)))
end

function matrixInfPNorm(A::Matrix{T}; p=2) where T <: Number
    jthBit = (i,j)->1&(i>>(j-1))
    enumeratedCols = collect(enumerate(eachcol(A)))
    colMultipleIter = map(
        i -> map(
            ((j,col),) -> (-2*jthBit(i,j)+1)*col,
            enumeratedCols
        ),
        0:2^(size(A,2)-1)-1
    )

    return maximum(map(cols -> norm(sum(cols), p), colMultipleIter))
end

#Polynomial Norms

function PolyNorm1(p::HCMK.Expression)
    return sum(abs.(HCMK.coefficients(p,HCMK.variables(p))))
end

function PolyNormW(p::HCMK.Expression)
    (M,C)=HCMK.exponents_coefficients(p,HCMK.variables(p))
    return sqrt(sum((C.^2)./Multinomials(M)))
end

function Multinomials(M::Matrix{S}) where S<:Integer
    return [ factorial(sum(a))/reduce(*,map(factorial,a)) for a in collect(eachcol(M))]
end

#GridPolynomialSystemConstructors

struct PolynomialSystem{T <: Number}
    evaluate::Function
    jacobian::Function
    degrees::Vector{Int}
    coefficients::Vector{Vector{T}} #to be eliminated
    monomialDegrees::Vector{Vector{Vector{Integer}}} #to be eliminated
end

function (p::PolynomialSystem{T})(x::Vector{T}) where T
    p.evaluate(x)
end

function GridSystem(polysys::HCMK.System,PolyNorm)

    invnorms = [
        1/PolyNorm(poly)
        for poly in polysys.expressions
    ]

    gridpolysys::PolynomialSystem{Float64} =
         PolynomialSystem{Float64}(
             v -> invnorms.*polysys(v),
             v -> Diagonal(invnorms).*HCMK.jacobian(polysys, v),
             HCMK.degrees(polysys),
             [],
             []
    )

    return gridpolysys
end

#function WeylRandomSystem

#function KacRandomSystem


#To delete...

function weylNorm(
    degree::Integer,
    monomials::Vector{Vector{S}},
    coefficients::Vector{T}
) where {S <: Integer, T <: Number}
    return sqrt(sum([
        _multinomial(degree, monomials[i])*coefficients[i]^2
        for i in eachindex(coefficients)
    ]))
end

function _multinomial(numerator::Integer, denominators::Vector{S}) where S<:Integer
    return factorial(numerator)/reduce(*,map(factorial,denominators))
end

function normalizePoly(P::PolynomialSystem{T}, polyNorm) where T
    invnorms = [
        1/polyNorm(P.degrees[i], P.monomialDegrees[i], P.coefficients[i])
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

end #module