## The utils.jl file contains utility functions for working with polynomials and their roots.

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
    # U, svals, Vt = LA.svd(M; alg = LA.QRIteration())
    svals = LA.svdvals(M) # Does not have option 'alg'.
    return last(svals)
    # return svals[findlast(>(ε), svals)]
end

normsqsq(A) = sum(a -> a^2, A)
Onorm(A::AbstractArray) = maximum(abs.(A))

homothecy(m) = x -> m.*x
translation(a) = x -> x.+a
dilation(a, m) = translation(a)∘homothecy(m)

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
    Id(::Type{T}, n)

Returns the identity matrix of size `n` and type `T`.
"""
function Id(::Type{T}, n) where T
  return LA.Diagonal(ones(T, n))
end

"""
    Id(iter)

Returns the identity matrix of size `length(iter)` and type `eltype(iter)`.
"""
Id(iter) = Id(eltype(iter), length(iter))
# diagonal(v) = cat(v...; dims = (1, 2)) # not using LinearAlgebra.jl

"""
    Δ(F)

Returns the degree matrix of a system of polynomials `F`.
"""
function Δ(F)
    return LA.Diagonal(degrees(F))
end

"""
    D(x)

Returns the matrix difference between the identity matrix and the outer product of `x` and its conjugate transpose.
"""
function D(x)
    # println(typeof(x), " ", x)
    return Id(x) - x*transpose(conj(x))
end
