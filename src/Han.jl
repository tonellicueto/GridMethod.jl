module Han
using LinearAlgebra
using ..GridModule
using ..GridBatch
using ..ConditionNumbers
using ..Norms
using .Iterators

export gridHan!
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
    while length(coordinates) > 0
        Threads.@threads for coord in coordinates 
            node = GridNode(G,depth,coord)
            if _HanCondition(node) || (!isnothing(maxDepth) && depth == maxDepth)
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