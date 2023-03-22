

divpow2(n, N) = n/exp2(N)
bound(depth) = 2^(depth)-1
function odds(depth::Int)
    bnd = bound(depth)
    return -bnd:2:bnd
end

function oddsgrid(dim, depth)
    return Iterators.product(ntuple(_ -> odds(depth), Val(dim))...)
end

function grid_getcoords(iter, depth)
    return Iterators.map(n -> divpow2(n, depth), iter)
end

# grid_getallcoords(depth) = grid_getcoords(odds(depth), depth)

# function get_gridnodes{dim}(depth)
#     return Iterators.map(I -> grid_getcoords(I, depth), oddsgrid(dim, depth))
# end

struct GridNode{T, dim}
    odds::Vector{Int}
    depth::Int
    image::T
end

image(g::GridNode) = g.image
function node(g::GridNode)
    @unpack odds, depth = g
    return grid_getcoords(odds, depth)
end

function point(g::GridNode, C = nothing)
    @unpack odds, depth, image = g
    if C == nothing
        return (odds, collect(grid_getcoords(odds, depth)), image)
    end
    return (odds, collect(grid_getcoords(odds, depth)), depth, C*inv(exp2(depth)), image)
end


mutable struct Grid{T, dim, F}
    nodes::Vector{GridNode{T, dim}}
    fun::F
end

nodes(G::Grid) = G.nodes
fun(G::Grid) = G.fun
dim(::Grid{T, d, F}) where {T, d, F} = d

Base.push!(G::Grid, g::GridNode) = push!(nodes(G), g)
Base.append!(G::Grid, iter) = append!(nodes(G), iter)
Base.pop!(G::Grid) = Base.pop!(nodes(G))
Base.isempty(G::Grid) = Base.isempty(nodes(G))
Base.findmin(G::Grid) = Base.findmin(image.(nodes(G)))

function points(G::Grid, C = nothing)
    @unpack nodes = G
    return [point(g, C) for g in nodes]
end

## Empty grid
Grid(::Type{T}, dim::Int, fun) where T = Grid{T, dim, typeof(fun)}(GridNode{T, dim}[], fun)
Grid(dim::Int, fun) = Grid(Float64, dim, fun)

## Generate nodes of a grid
GridNode(G::Grid, odds, depth, image) = GridNode{typeof(image), dim(G)}(collect(odds), Int(depth), image)

function GridNode(G::Grid, odds, depth)
    node = grid_getcoords(odds, depth)
    image = fun(G)(node...)
    GridNode(G, odds, depth, image)
end

function GridNodes(G::Grid, iter_odds, depth)
    Iterators.map(I -> GridNode(G, I, depth), iter_odds)
end

## Initial grid at some `depth`
function Grid(::Type{T}, dim::Int, fun, depth::Int) where T
    G = Grid(T, dim, fun)
    append!(G, GridNodes(G, oddsgrid(dim, depth), depth))
    return G
end
Grid(dim::Int, fun, depth::Int) = Grid(Float64, dim, fun, depth)

## Append in a grid `G` the next nodes (subdivisions) of a node `g`
function nextleaves(odds::Vector{Int})
    oodds = 2 .* odds
    col1 = oodds .+ 1
    col2 = oodds .- 1
    return Iterators.product(eachrow(hcat(col1, col2))...)
end

function append_nextleaves!(G::Grid, g::GridNode)
    @unpack odds, depth = g
    nextodds = nextleaves(odds)
    append!(G, GridNodes(G, nextodds, depth +1))
    return G
end

"""
    refine_grid(T=Float64, fun, C, m, dim; step₀ = Int(ceil(log2(C))),
                G = Grid(dim, fun, max(0, depth₀)), isfine = _isHan,
                append_subdivision! = append_nextleaves!)

Returns a list of `grid_point` objects that form a refined grid satisfying the Han condition.

# Arguments:
 - `fun`: The function to be evaluated at each point in the grid.
 - `C`: A parameter used in the Han condition.
 - `m`: Lower bound of the image.
 - `dim`: The dimension of the grid.

# Keyword arguments:
 - `depth₀`: Depth of the initial grid.
 - `G`: Initial grid, default is grid of depth `depth₀`.
 - `isfine`: A function that determines whether a point is valid in the grid.
 - `append_subdivision`: A function that appends to the grid the subdivision of a point.
"""
function refine_grid(fun, C, m, dim; depth₀ = Int(ceil(log2(C))),
                     G = Grid(dim, fun, max(0, depth₀)), isfine = _isHan,
                     append_subdivision! = append_nextleaves!)
    H = Grid(dim, fun)
    while !(isempty(G))
        g = pop!(G) # Removes an element of G and stores it at g.
        @unpack odds, depth, image = g
        if image < m
            @warn "Small norm" image, m
            return H # WARN Function returning two different types
        end
        if isfine(image, depth, C)
            # println("Step:", depth, " C:", C, "--", isfine(image, depth, C), image)
            push!(H, g)
        else
            append_subdivision!(G, g)
        end
    end
    return H
end


# All the evaluations of the function `fun` are done via this function.
# So, modify here for more efficient implementations.
# """
#     grid_point(fun::Function, x, step)

# Creates a new `grid_point` object with the given function, coordinates, and step, where the function value is computed by calling `fun(x...)`.
# """
# grid_point(fun::Function, x, step) = grid_point(x, fun(x...), step)

# function Grid{T, dim}(nodes, images, steps) where {T, dim}
#     U = promote_type(eltype(eltype(nodes)), eltype(images))

# end

# function root(G::Grid{T, d}, i) where {T, d}
#     new_data = TreeData{T, d}(nodes(G), step(G,i))
#     Tree{T, d}(; data = new_data )
# end

# step_to_factor(::Type{T}, step) where T = T(^(RATIO[], step)) # 1/2^i
# step_to_factor(step) = step_to_factor(Float64, step)

# """
#     step_to_ratio(step) = inv(exp2(step))

# Converts the subdivided index `step` of a `grid_point` into the ratio of such point.
# """
# step_to_ratio(step) = inv(exp2(step)) # 1/2^i

# """
#     radius(g::grid_point)

# Returns the radius of the `grid_point` object `g`, which is computed as `1/2^i`, where `i` is the step at which the point was added to the grid.
# """
# radius(g::grid_point) = step_to_ratio(g.step)

# """
#     eltype(::grid_point{T}) where {T}

# Returns the type of the elements in the `x` field of a `grid_point` object.
# """
# Base.eltype(::grid_point{T}) where T = T

# function Base.show(io::IO, point::grid_point)
#     print(io, "$(eltype(point))[")
#     print(io, join(point.x, ", "))
#     print(io, "], $(point.fx), $(point.step)")
# end



# function findallandnot(f::Function, A)
#     I = findall(f, A)
#     return I, setdiff(eachindex(A), I)
# end

# function push_subdivided!_or_fine!(H, G, br, isfine, m, C, fun)
#     g = pop!(H)
#     @unpack x, fx, step = g
#     if fx < m
#         @warn "Small norm" fx, m
#         return m # WARN Function returning two different types
#     end
#     if isfine(fx, step, C)
#         push!(G, g)
#     else
#         pushsubdivided!(H, fun, x, br)
#     end
# end

# # checknorms(g, m) = (g.fx < m)
# # checknorms(G, m) = any(broadcast(checknorms, G, m))

# function fine_grid(fun, C, m, dim; isfine = _isHan,
#                       pushsubdivide! = cube_pushsubdivided!)
#     I = map(isfine, G)
#     G, H = G[I], G[(!).(I)]
#     step = 1
#     while !(isempty(H))
#         if checknorms(H, m)
#             @warn "Small norm" fx, m
#             return m
#         end
#         br = branches(G, i) # Compute branches for the whole set of pt with step `i`.
#         subdivide!(H, G, m, br)
#     end
# end

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
#         # H = H[J]
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






#### Original code below:
#### Explanation and documentation generated with Chat GPT.
#### Steps:
####   1.- "Explain the following code written in Julia language used to compute
####        [explain something, use keywords]
####        ```julia
####        [code goes here]
####        ```"
####   2.- "Use the previous explanation to document de code.
####        Here it is an example of a documented function in Julia language:
####        ```julia
####        [example following your style here]
####        ```
####        "
####   3.- Copy-paste, remove extra sentences, rephrase some sentences... That's pretty much it!
# struct grid_point{T}
#     x::Vector{T}
#     fx::T
#     step::Int
# end
# function grid_point(x, fx, step)
#     U = promote_type(eltype(x), typeof(fx))
#     return grid_point{U}(convert(Vector{U}, x), convert(U, fx), Int(step))
# end
# grid_point(fun::Function, x, step) = grid_point(x, fun(x...), step)

# point(g::grid_point) = g.x
# fpoint(g::grid_point) = g.fx
# radius(g::grid_point) = step_to_ratio(g.step)
# Base.eltype(::grid_point{T}) where {T} = T

# function cube_nthinitgrid(fun, dim, n; ratio = 0.5, root = tree_root(dim),
#                           coord_opt = coord_opt(dim))
#     tree = tree_nthleaves(dim, n; ratio = ratio, root = root, coord_opt = coord_opt)
#     return [grid_point(fun, p, n) for p in tree]
# end

# function cube_pushsubdivided!(G, p, step, fun; dim = length(p),
#                               coord_opt = coord_opt(eltype(p), dim))
#     ratio = step_to_ratio(step + 1)
#     for q in tree_nextleaves(p, ratio; dim = dim, coord_opt = coord_opt)
#         push!(G, grid_point(fun, q, step + 1))
#     end
# end

# _isHan(fx, step, C) = C * step_to_ratio(step) < fx

# function refine_grid(fun, C, m, dim; step₀ = Int(ceil(log2(C))),
#                      G = cube_nthinitgrid(fun, dim, step₀), isfine = _isHan,
#                      pushsubdivided! = cube_pushsubdivided!)
#     println(step₀)
#     H = eltype(G)[]
#     while !(isempty(G)) # That is, while G is not the empty array.
#         g = pop!(G) # Removes an element of G and stores it at g.
#         @unpack x, fx, step = g
#         if fx < m
#             @warn "Small norm" fx, m
#             return m # WARN Function returning two different types
#         end
#         if isfine(fx, step, C)
#             push!(H, g)
#         else
#             pushsubdivided!(G, x, step, fun)
#         end
#     end
#     return H
# end
# function refine_grid!(G, fun, C, m, dim; isfine = _isHan,
#                       pushsubdivide! = cube_pushsubdivided!)
#     return G = refine_grid(fun, C, m, dim; step₀ = 0, G = G, isfine = isfine,
#                            pushsubdivide! = pushsubdivide!)
# end

# # Question: IMPLEMENT finding the min in H's 'while' building loop?
# function Han_min(H, L, C)
#     L < C || @warn "Required condition 0<1-inv(C)L not satisfied." L, C, inv(C) * L
#     γ = 1 - inv(C) * L # This is positive because precondition L < C
#     minfx, n = findmin(fpoint, H)
#     return γ * minfx, n
# end
