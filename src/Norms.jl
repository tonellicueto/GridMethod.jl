"""
    $(SIGNATURES)


    unfold_do_fold takes a polynomial system poly_sys and then
    1) unfolds the polynomials into terms
    2) applies the map term_map to each term_map
    3) puts together the values of the term of each polynomial using fold_poly (by default sum)
    4) puts together the values of the polynomials using fold_sys (by default sum)
"""

module Norms
using ..Polynomial
using .Iterators
using LinearAlgebra
using Logging

export polyNorm1
export hinfNorm
export matrixInfPNorm

function unfold_do_fold(poly_sys, term_map, fold_sys = sum, fold_poly = sum; kwargs...)
    poly_values = Iterators.map(poly_terms -> fold_poly(term_map, poly_terms),
                                 terms(poly_sys; kwargs...)) ## See ModelKitDocs.jl
    return fold_sys(poly_values)
end

"""
    Wnorm(F)
Weil norm of a system of polynomials.
"""

function Wnorm(F; kwargs...)
    function term_Wsqnorm(IC)
        ## Unpack exponent and coefficient
        (I, C) = IC
        return normsqsq(C)*inv(Combinatorics.multinomial(I...))
    end
    return sqrt(unfold_do_fold(F, term_Wsqnorm; kwargs...))
end

"""
    Onorm(F)

1-norm of a system of polynomials.
"""

function Onorm(F; kwargs...)
    return unfold_do_fold(F, abs∘last, maximum; kwargs...)
end

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
"""
THINGS TO DO:
1) Define normalization of a system:
Wnormalize(F) : divide each polynomial by its Wnorm
Onormalize(F) : divide each polynomial by its Onorm
2) Define
LInftynorm(F) : compute maximum of absolute value of each polynomial on sphere and take the maximum.
  -- heuristics: find something...
  -- construct grid and maximize over it.
and
LInftynormalize(F)
"""
