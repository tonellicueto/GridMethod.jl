import "Polynomial.jl"

function local_condition(F,x, norm_sys; degreematrix, get_sigma, norm_image, norm_denominator, kwargs...)
    norm₀ = norm_sys(F; kwargs...) #addway to input the precomputed norm
    x₀ = collect(promote(x...))
    (fx, Jfx) = eval_and_J(x₀, F; kwargs...) #add way to input the precomputed value and Jacobian
    σ = get_sigma(x₀, Jfx, degreematrix(F)) #add way to input the precomputed degreematrix
    normfx = norm_image(fx)
    return norm₀/norm_denominator(normfx, σ)
end

function local_condition_number(
    polysys::PolynomialSystem{T},
    x::Vector{T},
    
) where T
end

"""
    Id(::Type{T}, n)

Returns the identity matrix of size `n` and type `T`.

Id(iter)

Returns the identity matrix of size `length(iter)` and type `eltype(iter)`.
"""
function Id(::Type{T}, n) where T
  return LA.Diagonal(ones(T, n))
end

Id(iter) = Id(eltype(iter), length(iter))
# diagonal(v) = cat(v...; dims = (1, 2)) # not using LinearAlgebra.jl

"""
    tangent_proj(x)

Returns the matrix difference between the identity matrix and the outer product of `x` and its conjugate transpose.
"""
function tangent_proj(x)
    # println(typeof(x), " ", x)
    return Id(x) - x*transpose(conj(x))
end

"""
   sqDelta(F)

Returns the degree matrix of a system of polynomials `F`.
"""
function sqDelta(F)
    return LA.sqrt(LA.Diagonal(degrees(F))) #simplify this...
end

"""
    K(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""

Wcond(F,x; kwargs...) = local_condition(F,x, Wnorm;
                              degreematrix = sqDelta,
                              get_sigma = lastsingvaluesq,
                              norm_image = normsqsq,
                              norm_denominator = (x,y)->sqrt(x+y),
                              kwargs...)

lastsingvaluesq(x, Jfx, sqΔ) = last_sval(inv(sqΔ)*Jfx*tangent_proj(x))^2

# function K(Wnorm, Δ, x, fx, Jx)
#     σ = get_sigmaK(x, Jx, Δ)
#     f2 = normsqsq(fx)
#     return Wnorm/sqrt(f2 + σ)
# end


"""
    Delta(F)

Returns the degree matrix of a system of polynomials `F`.
"""
function Delta(F)
    return LA.Diagonal(degrees(F))
end

"""
    C(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
Ocond(F, x; kwargs...) = local_condition(F, x, Onorm;
                              degreematrix = Delta,
                              get_sigma = lastsigma_inf,
                              norm_image = fx -> maximum(abs.(fx)),
                              norm_denominator = max,
                              kwargs...)

lastsigma_inf(x, Jfx, Δ) = 1/LA.opnorm(inv(Jfx)*Δ,Inf)

# function C(Onorm, Δ)
#     (xfJ) -> begin
#         # unpack
#         (x, fx, Jx) = xfJ
#         # M = inv(Jx*D(x))*Δ
#         M = (Jx*D(x))*inv(Δ)
#         # σ = inv(maximum(M))
#         σ = maximum(M)
#         f2 = maximum(abs.(fx))
#         return Onorm/max(f2, σ)
#     end
# end


