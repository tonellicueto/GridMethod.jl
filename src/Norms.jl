module Norms
using ..Polynomial
using .Iterators
using LinearAlgebra
using Logging

export polyNorm1
export hinfNorm
export matrixInfPNorm

function polyNorm1(poly::PolynomialSystem{T})::T where T
    # the 1-norm of a polynomial system is the maximum of
    # the absolute sum of the coefficients from every polynomial
    # in the system.
    return maximum(map(sum, map(l -> map(abs, l),poly.coefficients)))
end

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
end # module