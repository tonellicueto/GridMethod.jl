## The utils.jl file contains utility functions for working with polynomials.

"""
    rand_poly(::Type{T}=Float64, X, d; coeffs = I -> randn(T), n = length(X))

Returns a random polynomial of coefficient-type `T` in variables `X` with degree `d`.

# Keyword arguments:
 - `coeffs`: A function that takes a tuple of exponents and returns a coefficient.
 - `n`: The number of variables.
"""
function rand_poly(::Type{T}, X, d; coeffs = I -> randn(T), n = length(X), homo = false) where T
    if homo
        expos = Combinatorics.multiexponents(n, d)
    else
        expos = Iterators.map(x -> x[begin+1:end], Combinatorics.multiexponents(n+1, d))
    end
    sum(I -> coeffs(I)*monomial(X, I), expos)
end
rand_poly(X, d; coeffs = I -> randn(Float64), n = length(X), homo = false) =
    rand_poly(Float64, X, d; coeffs = coeffs, n = n, homo = homo)

const ε = 1e-6 # bound for meaning zero
function last_sval(M)
    U, svals, Vt = LA.svd(M; alg = LA.QRIteration())
    # svals = LA.svdvals(M) # Does not have the keyword 'alg'.
    return last(svals)
    # return svals[findlast(>(ε), svals)]
end

normsqsq(A) = sum(a -> a^2, A)
norm1(A::AbstractArray) = maximum(abs.(A))

# pmones(::Type{T}, ::Val{d}) where {T, d} =
#     Iterators.map(collect,
#                   Iterators.ProductIterator(
#                       ntuple(_ -> (one(T), -one(T)),
#                              Val(d),
#                              )))

homothecy(m) = x -> m.*x
translation(a) = x -> x.+a
dilation(a, m) = translation(a)∘homothecy(m)

function terms(sys;
               vars = get_vars(sys),
               expanded = true)
    return  terms.(get_polys(sys);
                   vars = vars,
                   expanded = expanded)
end

## From HC.jl
"""
    monomial(X, I, n = length(X))

Returns the monomial of variables `X` with exponents `I`.

# Keyword arguments:
 - `n`: The number of variables.
"""
function monomial(X, I, n=length(X))
    prod(i -> X[i]^I[i], 1:n)
end

"""
    Id(::Type{T}=Float64, n)

Returns the identity matrix of size `n` and type `T`.
"""
function Id(::Type{T}, n) where T
    LA.Diagonal(ones(T, n))
end
Id(v::AbstractVector) = Id(eltype(v), length(v))
# diagonal(v) = cat(v...; dims = (1, 2)) # Without LinearAlgebra.jl

"""
    Δ(F)

Returns the degree matrix of a system of polynomials `F`.
"""
function Δ(F)
    LA.Diagonal(degrees(F))
end

"""
    D(x)

Returns the matrix difference between the identity matrix and the outer product of `x` and its conjugate transpose.

# Keyword arguments:
 - `n`: The size of the identity matrix.
"""
function D(x)
    Id(x) - x*transpose(conj(x))
end

