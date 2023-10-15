function local_condition(F,x, norm_sys; degreematrix, get_sigma, norm_image, norm_denominator, kwargs...)
    norm₀ = norm_sys(F; kwargs...)
    x₀ = collect(promote(x...))
    (fx, Jfx) = eval_and_J(x₀, F; kwargs...)
    σ = get_sigma(x₀, Jfx, degreematrix(F))
    normfx = norm_image(fx)
    return norm₀/norm_denominator(normfx, σ)
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
    Δ(F)

Returns the degree matrix of a system of polynomials `F`.
"""
function sqDelta(F)
    return sqrt(LA.Diagonal(degrees(F)))
end

"""
    K(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""

local_Wcond(F,x; kwargs...) = condition(F,x, Wnorm;
                              degreematrix = sqDelta,
                              get_sigma = lastsingvaluesq,
                              norm_image = normsqsq,
                              choosing = (x,y)->sqrt(x+y),
                              kwargs...)

lastsingvaluesq(x, Jfx, sqΔ) = last_sval(inv(sqΔ)*Jfx*tangent_proj(x))^2

# function K(Wnorm, Δ, x, fx, Jx)
#     σ = get_sigmaK(x, Jx, Δ)
#     f2 = normsqsq(fx)
#     return Wnorm/sqrt(f2 + σ)
# end


"""
    Δ(F)

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
Ocond(x, F; kwargs...) = condition(x, F, Onorm;
                              degreematrix = Delta,
                              get_sigma = lastsigma_inf,
                              norm_image = fx -> maximum(abs.(fx)),
                              choosing = max,
                              kwargs...)

lastsigma_inf(x, Jfx, Δ) = 1/maximum(Onorm, eachrow(inv(Δ)*inv(Jfx)))

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
