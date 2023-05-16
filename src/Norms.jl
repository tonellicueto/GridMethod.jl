
"""
    $(SIGNATURES)

Returns the norm of the system of polynomials `F` given by applying `f` to terms `IC = (I, C)` and folding each polynomial accordingly to `fold_poly` and finally folding the system by `fold_system`.
"""
function norm(sys, fun_term, fold_system = sum, fold_poly = sum; kwargs...)
    poly_norms = Iterators.map(poly_terms -> fold_poly(fun_term, poly_terms),
                                 terms(sys; kwargs...)) ## See ModelKitDocs.jl
    return fold_system(poly_norms)
end

function Wnorm_f(IC)
    ## Unpack exponent and coefficient
    (I, C) = IC
    return normsqsq(C)*inv(Combinatorics.multinomial(I...))
end

"""
    Wnorm(F)
Weil norm of a system of polynomials.
"""
function Wnorm(F; kwargs...)
    return sqrt(norm(F, Wnorm_f; kwargs...))
end

"""
    norm1(F)

L_1 norm of a system of polynomials.
"""
function norm1(F; kwargs...)
    return norm(F, abs∘last, maximum; kwargs...)
end
