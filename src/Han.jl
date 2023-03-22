

"""
    _isHan(image, depth, C)

Returns the boolean value of `C * inv(exp2(depth)) < image`.
"""
_isHan(image, depth, C) = C * inv(exp2(depth)) < image
# _isHan(fx, step, C) = C * step_to_ratio(step) < fx

"""
    Han_min(H, L, C)

Returns the minimum function value in `H` scaled by the factor `γ = 1 - inv(C) * L`.
"""
function Han_min(G, L, C)
    L < C || @warn "Required condition 0<1-inv(C)L not satisfied." L, C, inv(C) * L
    minfx, n = findmin(G)
    # \gamma is positive by precondition L < C
    return (1 - inv(C) * L) * minfx, n
end

