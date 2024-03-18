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
            if _HanCondition(node)
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

function _HanCondition(node::GridNode{T, dim}) where {T, dim}
    return condition(node)/2^depth(node)≤0.5
end

function _splitNode(node::GridNode{T, dim}) where {T, dim}
    depth = node.depth + 1
    center = coordinates(node)
    return map(
        t -> [n/2^depth for n in t] + center,
        product(repeated([-one(T),one(T)],dim)...)
    )
end
end #module