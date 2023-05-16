
"""
    grid_Han(T=Float64, fun, C, m, dim;
                depth₀ = Int(ceil(log2(C))))

Returns a grid satisfying the Han condition.

# Arguments:
 - `T`: eltype of grid nodes
 - `fun`: The function to be evaluated at each point in the grid.
 - `C`: A parameter used in the Han condition.
 - `m`: Lower bound of the image.
 - `dim`: The dimension of the grid.

# Keyword arguments:
 - `depth₀`: Depth of the initial grid.
 - `G₀`: Initial grid, default is grid of depth `depth₀`.
"""
function grid_Han(::Type{T}, fun, C, m, dim;
                     depth₀ = Int(ceil(log2(C)))) where T
    H = Grid(T, dim, fun) ## Empty grid
    G = Grid(T, dim, fun, max(0, depth₀)) ## Grid of min depth.
    while !(isempty(G))
        g = pop!(G) # Removes an element of G and stores it at g.
        @unpack odds, depth, image = g
        if image < m
            @warn "Small norm" image, m
            return H # WARN Function returning two different types
        end
        if C/exp2(depth) < image
            # println("Step:", depth, " C:", C, "--", isHan(image, depth, C), image)
            push!(H, g)
        else
            append_nextleaves!(G, g)
        end
    end
    return H
end

grid_Han(fun, C, m, dim;
         depth₀ = Int(ceil(log2(C)))) = grid_Han(Float64, fun, C, m, dim;
                                                 depth₀ = depth₀)

Han_factor(C, L) = 1 - L/C

"""
    Han_min(G, L, C)

Returns the minimum function value in grid `G` scaled by the factor `γ = 1 - L/C`.
"""
function Han_findmin(G, L, C)
    L < C || @warn "Required condition 0<1-L/C not satisfied." L, C, L/C
    minfx, n = findmin(G)
    print("γ = $(Han_factor(C, L))\nimage = $(image(minfx))\n")
    return Han_factor(C, L)*image(minfx), n
end

# function refine_grid!(G, fun, C, m, dim; isfine = _isHan,
#                       pushsubdivide! = cube_pushsubdivided!)
#     if checknorms(G, m)
#         @warn "Small norm" fx, m
#         return m
#     end
#     I = map(isfine, G)
#     G, H = G[I], G[(!).(I)]
#     while !(isempty(H))
#         br = branches(G, H)
#         # subdivide!(H, br)
#         # if checknorms(H, m)
#         #     @warn "Small norm" fx, m
#         #     return m
#         # end
#         # I = map(isfine, H)
#         # append!(G, H[I])
#         # H = H[(!).(I)]
#         for _ in 1:length(H)
#             g = popfirst!(H)
#             @unpack x, fx, step = g
#             if fx < m
#                 @warn "Small norm" fx, m
#                 return m
#             end
#             if isfine(fx, step, C)
#                 push!(G, g)
#             else
#                 pushsubdivided!(H, fun, x, br[step])
#             end
#         end
#     end
# end

# eachstep(H) = unique(map(h -> h.step, H))
# branches(G, H) = Dict([(i, branches(G,i)) for i in eachstep(H)])
