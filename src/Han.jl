module Han
using LinearAlgebra
using ..GridModule
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
    maxDepth::Union{UInt,Nothing}=nothing,
    offset::Union{Vector{T},Nothing}=nothing,
    _split=splitCoordinate
) where {T, dim}
    local_depth = depth
    coordinates = collect(generateCoordinates(
        T,
        dim,
        local_depth,
        lower,
        upper
    ))
    if !isnothing(offset)
        broadcast!(v -> v+offset,coordinates,coordinates)
    end
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
                newCoordinates = _split(
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

function projectiveGridHan!(
    PG::ProjectiveGrid{T, dim},
    depth::UInt;
    lower::T=-1*one(T),
    upper::T=one(T),
    maxDepth::Union{UInt,Nothing}=nothing
) where {T, dim}
    basisVectors = usualBasis(T, UInt(dim))
    Threads.@threads for (i,G) in enumerate(PG.charts)
        basisIndex = ((i-1)%dim)+1
        vectorSign = (-1)^(floor(Int,(i-1)/dim))
        gridHan!(
            G,
            depth;
            lower=lower,
            upper=upper,
            maxDepth=maxDepth,
            offset=basisVectors[basisIndex]*vectorSign,
            _split=_projectiveSplitCoordinate(basisIndex)
        )
    end
    
    PG.est_condition = maximum(est_condition, PG.charts)
end

function _projectiveSplitCoordinate(dimExclude::UInt)
    function __projectiveSplitCoordinate(
        v::Vector{T},
        target_depth::UInt,
        depth::UInt,
        scale::T
    )
        vProjection = cat(v[1:dimExclude-1],v[dimExclude+1:length(v)])
        coords = splitCoordinate(
            vProjection,
            target_depth;
            depth=depth,
            scale=scale
        )
        return [
            cat(u[1:dimExclude-1],v[dimExclude],u[dimExclude:length(u)])
            for u in coords
        ]
    end

    return __projectiveSplitCoordinate
end
end #module
