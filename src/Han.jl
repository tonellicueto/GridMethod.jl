module Han
using LinearAlgebra
using ..GridModule
using ..GridBatch
using ..ConditionNumbers
using ..Norms
using .Iterators

export gridHan!
##TODO: Rewrite Han to directly compute the grid.
## It should initiate computing the degree matrix, the norm of each polynomial in the system and
## normalize each polynomial in the system by its norm.
## It should start subdividing storing the value, jacobian, and condition at
## each point of the grid.
## It should output the final grid, together with the upper estimate of the condition number and the normalized system.
## One can write a subprogram to normalize the system, i.e., divide each polynomial by its norm.
function gridHan!(
    G::Grid{T, dim},
    depth::UInt;
    lower::T=-1*one(T),
    upper::T=one(T),
    maxDepth::Union{UInt,Nothing}=nothing
) where {T, dim}
    coordinates = collect(generateCoordinates(
        T,
        dim,
        depth,
        lower,
        upper
    ))
    reprocess::Vector{Vector{T}} = []
    reprocessLock = ReentrantLock()
    gridLock = ReentrantLock()
    while length(coordinates) > 0 && (isnothing(maxDepth) || depth <= maxDepth)
        Threads.@threads for coord in coordinates 
            node = GridNode(G,depth,coord)
            if norm(image(node),Inf)*condition(node)≥1
                lock(gridLock)
                try
                    push!(G,node)
                finally
                    unlock(gridLock)
                end
            else
                newCoordinates = _splitNode(node)
                lock(reprocessLock)
                try
                    append!(reprocess, newCoordinates)
                finally
                    unlock(reprocessLock)
                end
            end
        end

        coordinates = reprocess
        reprocess = []
        depth = depth + 1
    end
end

function _splitNode(node::GridNode{T, dim}) where {T, dim}
    depth = node.depth + 1
    center = coordinates(node)
    return map(
        t -> [n/2^depth for n in t] + center,
        product(repeated([-one(T),one(T)],dim)...)
    )
end

#function grid_Han(::Type{T}, objective_to_minimize, C, m, dim;
#                     depth₀ = Int(ceil(log2(C)))) where T
#    H = Grid(T, dim, objective_to_minimize) ## Empty grid
#    G = Grid(T, dim, objective_to_minimize, max(0, depth₀)) ## Grid of min depth.
#    while !(isempty(G))
#        g = pop!(G) # Removes an element of G and stores it at g.
#        @unpack odds, depth, image = g
#        if image < m
#            @warn "Small norm" image, m
#            return H # WARN Function returning two different types
#        end
#        if C/exp2(depth) < image
#            # println("Step:", depth, " C:", C, "--", isHan(image, depth, C), image)
#            push!(H, g)
#        else
#            append_nextleaves!(G, g)
#        end
#    end
#    return H
#end
#
#grid_Han(objective_to_minimize, C, m, dim;
#         depth₀ = Int(ceil(log2(C)))) = grid_Han(Float64, objective_to_minimize, C, m, dim;
#                                                 depth₀ = depth₀)
#
#Han_factor(C, L) = 1 - L/C
#
"""
    Han_min(G, L, C)

Returns the minimum function value in grid `G` scaled by the factor `γ = 1 - L/C`.
"""
#function Han_findmin(G, L, C)
#    L < C || @warn "Required condition 0<1-L/C not satisfied." L, C, L/C
#    minfx, n = findmin(G)
#    print("γ = $(Han_factor(C, L))\nimage = $(image(minfx))\n")
#    return Han_factor(C, L)*image(minfx), n
#end

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
end #module