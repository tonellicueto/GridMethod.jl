
## Utils for Grid
function grid_oddcoords(depth::Int)
    bound = 2^(depth)-1
    return -bound:2:bound
end

function grid_oddpoints(::Val{dim}, depth) where dim
    return Iterators.product(ntuple(_ -> grid_oddcoords(depth), Val(dim))...)
end

function grid_getcoords(depth, iter)
    return map(n -> n/exp2(depth), iter)
end

## A grid node
struct GridNode{T, dim}
    depth::Int
    odds::Vector{Int}
    image::T
end

depth(g::GridNode) = g.depth
image(g::GridNode) = g.image
dim(::GridNode{T, dim}) where {T, dim} = dim
function node(g::GridNode)
    @unpack depth, odds = g
    return grid_getcoords(depth, odds)
end

#### Grid interface
mutable struct Grid{T, dim, F}
    gridnodes::Vector{GridNode{T, dim}}
    fun::F
end

## Get float coordinates
nodes(G::Grid) = map(node, G)

## Access to grid
fun(G::Grid) = G.fun
gridnodes(G::Grid) = G.gridnodes
dim(::Grid{T, dim, F}) where {T, dim, F} = dim

## We are allowing the basic functions of julio to work in new data types
Base.length(G::Grid) = Base.length(gridnodes(G))
Base.size(G::Grid) = (Base.length(G),)
Base.IndexStyle(::Type{<:Grid}) = IndexLinear()
Base.getindex(G::Grid, i::Int) = getindex(gridnodes(G), i)


# ## Mutate a grid
Base.push!(G::Grid, g::GridNode) = Base.push!(gridnodes(G), g)
Base.pop!(G::Grid) = Base.pop!(gridnodes(G))

# ## Iterate over a grid
Base.eltype(G::Grid) = Base.eltype(gridnodes(G))
Base.isempty(G::Grid) = Base.isempty(gridnodes(G))
function Base.iterate(G::Grid, state = 1)
    return Base.iterate(gridnodes(G), state)
end

## So, we can call `maximum(G)` for a grid `G`
Base.isless(g::GridNode, h::GridNode) = Base.isless(image(g), image(h))

## So we can call `findmin(G)` for a grid `G`
Base.keys(G::Grid) = Base.keys(gridnodes(G))
Base.values(G::Grid) = Base.values(gridnodes(G))

function Base.print(io::IO, g::GridNode, kwargs...)
    @unpack depth, odds, image = g
    # print(io, "GridNode: \n")
    print(io, "┌Node: ", odds, " ", node(g), "\n")
    print(io, "└Depth: ", depth, " Image: ", image)
    print(io, kwargs...)
end

Base.show(io::IO, ::MIME"text/plain", g::GridNode) = Base.print(io, g)

function Base.show(io::IO, ::MIME"text/plain", G::Grid)
    for g in G
        print(io, g, "\n")
    end
end


## Given a grid `G`, create new grid nodes for `G`
# Take care of types
function GridNode(::Grid{T, dim, F}, depth, odds, image) where {T, dim, F}
    GridNode{T, dim}(Int(depth), collect(odds), T(image))
end

# Compute the image
function GridNode(G::Grid, depth, odds)
    float_coord = grid_getcoords(depth, odds)
    image = fun(G)(float_coord...)
    GridNode(G, depth, odds, image)
end

## Create an empty grid
function Grid(::Type{T}, dim::Int, fun) where T
    F = typeof(fun)
    gridnodes = GridNode{T, dim}[]
    return Grid{T, dim, F}(gridnodes, fun)
end
Grid(dim::Int, fun) = Grid(Float64, dim, fun)
# Base.empty!(G::Grid{T, dim, F}) where {T, dim, F} = gridnodes(G) = GridNode{T, dim}[]
Base.empty(G::Grid{T, dim, F}) where {T, dim, F} = Grid(T, dim, fun(G))

## Create full grid of depth `depth`
function Grid(::Type{T}, dim::Int, fun, depth::Int) where T
    G = Grid(T, dim, fun) ## Empty grid
    for oddpoint in grid_oddpoints(Val(dim), depth)
        push!(G, GridNode(G, depth, oddpoint))
    end
    return G
end
Grid(dim::Int, fun, depth::Int) = Grid(Float64, dim, fun, depth)

## Append in a grid `G` the next nodes (subdivisions) of a node `g`
# Given the vectors of odds `odds` representing `g`, returns
# an iterator over the divisions of `g` into `2^dim` (where `dim == length(odds) = true`)
# new points, one for each orthant.
function nextleaves(odds::Vector)
    oodds = 2 .* odds
    col1 = oodds .+ 1
    col2 = oodds .- 1
    return Iterators.product(eachrow(hcat(col1, col2))...)
end

function append_nextleaves!(G::Grid, g::GridNode)
    @unpack depth, odds = g
    for oddpoint in nextleaves(odds)
        push!(G, GridNode(G, depth +1, oddpoint))
    end
    return G
end


function grid_groupbydepth(G::Grid{T, dim, F}) where {T, dim, F}
    I = unique(depth.(G))
    out = Dict(I .=> [empty(G) for _ in I ])
    for g in G
        d = depth(g)
        push!(out[d], g)
    end
    return out
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
