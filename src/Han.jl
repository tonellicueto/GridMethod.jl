module Han
using LinearAlgebra
using ..GridModule
using ..GridBatch
using ..ConditionNumbers
using ..Norms
using ..Coordinates
using .Iterators

export gridHan!
function gridHan!(
    G::Grid{T, dim},
    depth::UInt;
    lower::T=-1*one(T),
    upper::T=one(T),
    maxDepth::Union{UInt,Nothing}=nothing
) where {T, dim}
    local_depth = depth
    coordinates = collect(generateCoordinates(
        T,
        dim,
        local_depth,
        lower,
        upper
    ))
    reprocess::Vector{Vector{T}} = []
    reprocessLock = ReentrantLock()
    gridLock = ReentrantLock()
    while length(coordinates) > 0
        Threads.@threads for coord in coordinates 
            node = GridNode(G,local_depth,coord)
            if _HanCondition(node) || (!isnothing(maxDepth) && local_depth == maxDepth)
                lock(gridLock)
                try
                    push!(G,node)
                finally
                    unlock(gridLock)
                end
            else
                newCoordinates = splitCoordinate(
                    GridModule.coordinates(node),
                    local_depth + 1;
                    depth=local_depth,
                    scale=upper-lower
                )
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
        local_depth += 1
    end

    # fill out estimated condition number
    G.est_condition = maximum(condition, G)
end

function _HanCondition(node::GridNode{T, dim}) where {T, dim}
    return condition(node)/2^depth(node)≤0.5
end
end #module