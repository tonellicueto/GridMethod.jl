
function constant(x, F, norm_sys; get_sigma, norm_image, choosing, kwargs...)
    norm₀, Δ₀ = norm_sys(F; kwargs...), Δ(F)
    (fx, Jx) = eval_and_J(collect(x), F; kwargs...)
    σ = get_sigma(x, Jx, Δ₀)
    f2 = norm_image(fx)
    return norm₀/choosing(f2, σ)
end

"""
    K(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
K(x, F; kwargs...) = constant(x, F, Wnorm;
                              get_sigma = get_sigmaK,
                              norm_image = normsqsq,
                              choosing = (x,y)->sqrt(x+y),
                              kwargs...)

get_sigmaK(x, Jx, Δ) = last_sval(inv(Δ)*Jx*D(x))

# function K(Wnorm, Δ, x, fx, Jx)
#     σ = get_sigmaK(x, Jx, Δ)
#     f2 = normsqsq(fx)
#     return Wnorm/sqrt(f2 + σ)
# end


"""
    C(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
C(x, F; kwargs...) = constant(x, F, norm1;
                              get_sigma = get_sigmaC,
                              norm_image = fx -> maximum(abs.(fx)),
                              choosing = max,
                              kwargs...)

get_sigmaC(x, Jx, Δ) = maximum(norm1, eachcol(inv(Δ)*Jx*D(x)))

# function C(norm1, Δ)
#     (xfJ) -> begin
#         # unpack
#         (x, fx, Jx) = xfJ
#         # M = inv(Jx*D(x))*Δ
#         M = (Jx*D(x))*inv(Δ)
#         # σ = inv(maximum(M))
#         σ = maximum(M)
#         f2 = maximum(abs.(fx))
#         return norm1/max(f2, σ)
#     end
# end
