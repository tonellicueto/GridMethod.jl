
"""
    norm(F, f, fold_system = sum, fold_poly = sum)

Returns the norm of the system of polynomials `F` given by applying `f` to terms `(I, C)` and folding each polynomial accordingly to `fold_poly` and finally folding the system by `fold_system`.
Assumes polynomials in `F` are expanded.
"""
function norm(F, f, fold_system = sum, fold_poly = sum; kwargs...)
    return fold_system(terms -> fold_poly(f, terms), terms(F; kwargs...)) ## trems in Utils.jl
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
    return norm(F, Wnorm_f; kwargs...)
end

"""
    norm1(F)

L_1 norm of a system of polynomials.
"""
function norm1(F; kwargs...)
    return norm(F, abs∘last, maximum; kwargs...)
end
