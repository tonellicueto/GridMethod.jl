module GridModule
using Parameters
using ..Polynomial
using ..ConditionNumbers
import Base: ==

export GridNode
export Grid
export depth
export coordinates
export image
export jacobian
export condition
export dim
export polysys
export gridnodes
export est_condition
export dim

struct GridNode{T <: Number, dim}
    depth::UInt
    coordinates::Vector{T}
    image::Vector{T}
    jacobian::Union{Matrix{T}, Nothing}
    condition::Union{T, Nothing}
end

#Functions for GridNode
depth(g::GridNode) = g.depth
coordinates(g::GridNode) = g.coordinates
image(g::GridNode) = g.image
jacobian(g::GridNode) = g.jacobian
condition(g::GridNode) = g.condition
dim(g::GridNode) = typeof(g).parameters[2]

function ==(node1::GridNode{T, dim}, node2::GridNode{T, dim}) where {T, dim}
    return (
        node1.depth == node2.depth
        && node1.coordinates == node2.coordinates
        && node1.image == node2.image
        && node1.jacobian == node2.jacobian
        && node1.condition == node2.condition
    )
end

# Structure defining the grid
mutable struct Grid{T <: Number, dim}
    polysys::PolynomialSystem{T}
    gridnodes::Vector{GridNode{T, dim}}
    est_condition::Union{T, Nothing}
end

#Functions to access Grid
polysys(G::Grid) = G.polysys
gridnodes(G::Grid) = G.gridnodes
est_condition(G::Grid) = G.est_condition
dim(G::Grid) = typeof(G).parameters[2] 

# Extension of Base functions of julia to Grid Data Types
Base.length(G::Grid) = Base.length(gridnodes(G))
Base.size(G::Grid) = (Base.length(G),)
Base.IndexStyle(::Type{<:Grid}) = IndexLinear()
Base.getindex(G::Grid, i::Int) = getindex(gridnodes(G), i)

# Base function to put a node in and a node out
Base.push!(G::Grid, g::GridNode) = Base.push!(gridnodes(G), g) #Add node g to grid G
Base.pushfirst!(G::Grid, g::GridNode) = Base.pushfirst!(gridnodes(G), g) #Add node g to grid G
Base.pop!(G::Grid) = Base.pop!(gridnodes(G)) #Picks a node g from G and removes it
Base.popfirst!(G::Grid, g::GridNode) = Base.popfirst!(gridnodes(G)) #Add node g to grid G

# Iterate over a grid
Base.eltype(G::Grid) = Base.eltype(gridnodes(G))
Base.isempty(G::Grid) = Base.isempty(gridnodes(G))
function Base.iterate(G::Grid, state = 1)
    return Base.iterate(gridnodes(G), state)
end

# So, we can call `maximum(G)` for a grid `G` to compute G.est_condition
Base.isless(g::GridNode, h::GridNode) = Base.isless(condition(g), condition(h))

# So we can call `findmax(G)` for a grid `G`
Base.keys(G::Grid) = Base.keys(gridnodes(G))
Base.values(G::Grid) = Base.values(gridnodes(G))

function Base.empty(G::Grid{T, dim}) where {T, dim}
    return Grid{T, dim}(polysys(G), [], est_condition(G))
end

function GridNode(
    grid::Grid{T, dim},
    depth::UInt,
    coordinates::Vector{T}
) where {T, dim}
    nodeImage = polysys(grid)(coordinates)
    nodeJacobian = polysys(grid).jacobian(coordinates)
    GridNode{T, dim}(
        depth,
        coordinates,
        nodeImage,
        nodeJacobian,
        localC(
            polysys(grid),
            coordinates;
            image=nodeImage,
            jacobian=nodeJacobian
        )
    )
end
end