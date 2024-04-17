module Han
using LinearAlgebra
using ..GridModule
using ..ConditionNumbers
using ..NormsPolynomials
using ..Coordinates
using .Iterators

export gridHan!
export projectiveGridHan!
export increaseMinDepth!
export increaseDepth!
export LazyIncreaseDepth!
export projectiveIncreaseMinDepth!
export projectiveIncreaseDepth!

function projectiveIncreaseMinDepth!(
    PG::ProjectiveGrid{T, dim},
    targetDepth::UInt;
    scale::T = one(T),
    nodeCondition=projectiveLocalC,
    nodeFilter=(node)->true
) where {T, dim}
    enumeratedCharts = collect(enumerate(PG.charts))
    Threads.@threads for (i,G) in enumeratedCharts
        basisIndex = UInt(((i-1)%dim)+1)
        increaseMinDepth!(
            G,
            targetDepth;
            scale=scale,
            nodeCondition=nodeCondition,
            nodeFilter=nodeFilter,
            _split=_projectiveSplitCoordinate(T, basisIndex)
        )
    end
end

function projectiveIncreaseDepth!(
    PG::ProjectiveGrid{T, dim},
    extraDepth::UInt;
    scale::T = one(T),
    nodeCondition=projectiveLocalC,
    nodeFilter=(node)->true
) where {T, dim}
    enumeratedCharts = collect(enumerate(PG.charts))
    Threads.@threads for (i,G) in enumeratedCharts
        basisIndex = UInt(((i-1)%dim)+1)
        increaseDepth!(
            G,
            extraDepth;
            scale=scale,
            nodeCondition=nodeCondition,
            nodeFilter=nodeFilter,
            _split=_projectiveSplitCoordinate(T, basisIndex)
        )
    end
end

function increaseMinDepth!(
    G::Grid{T, dim},
    targetDepth::UInt;
    scale::T = one(T),
    nodeCondition=localC,
    nodeFilter=(node)->true,
    _split=splitCoordinate
) where {T, dim}
    gridLock = ReentrantLock()
    Threads.@threads for _ in 1:length(G)
        node = nothing
        lock(gridLock)
        try
            node = popfirst!(G)
        finally
            unlock(gridLock)
        end
        newNodes = []
        if nodeFilter(node)
            if depth(node) ≥ targetDepth
                push!(newNodes,node)
            else
                newCoordinates = _split(
                    coordinates(node),
                    targetDepth;
                    depth=depth(node),
                    scale=scale
                )
                for coord in newCoordinates
                    newNode = GridNode(
                        G,
                        targetDepth,
                        coord;
                        localCondition=nodeCondition
                    )
                    if nodeFilter(newNode)
                        push!(newNodes,newNode)
                    end
                end
            end
        end 
        lock(gridLock)
        try
            append!(G,newNodes)
        finally
            unlock(gridLock)
        end
    end
end

function increaseDepth!(
    G::Grid{T, dim},
    extraDepth::UInt;
    scale::T = one(T),
    nodeCondition=localC,
    nodeFilter=(node)->true,
    _split=splitCoordinate
) where {T, dim}
    gridLock = ReentrantLock()
    Threads.@threads for _ in 1:length(G)
        node = nothing
        lock(gridLock)
        try
            node = popfirst!(G)
        finally
            unlock(gridLock)
        end
        newNodes = []
        if nodeFilter(node)
            newDepth = depth(node)+extraDepth
            newCoordinates = _split(
                coordinates(node),
                newDepth;
                depth=depth(node),
                scale=scale
            )
            for coord in newCoordinates
                newNode = GridNode(
                    G,
                    newDepth,
                    coord;
                    localCondition=nodeCondition
                )
                if nodeFilter(newNode)
                    push!(newNodes,newNode)
                end
            end
        end
        lock(gridLock)
        try
            append!(G,newNodes)
        finally
            unlock(gridLock)
        end
    end
end

function LazyIncreaseDepth!(
    G::Grid{T, dim},
    extraDepth::UInt;
    scale::T = one(T),
    nodeFilter=(node)->true,
    _split=splitCoordinate
) where {T, dim}
    gridLock = ReentrantLock()
    Threads.@threads for _ in 1:length(G)
        node = nothing
        lock(gridLock)
        try
            node = popfirst!(G)
        finally
            unlock(gridLock)
        end
        newNodes = []
        if nodeFilter(node)
            newDepth = depth(node)+extraDepth
            cond=condition(node)
            newCoordinates = _split(
                coordinates(node),
                newDepth;
                depth=depth(node),
                scale=scale,
            )
            for coord in newCoordinates
                newNode = LazyGridNode(
                    G,
                    newDepth,
                    coord,
                    cond
                )
                if nodeFilter(newNode)
                    push!(newNodes,newNode)
                end
            end
        end
        lock(gridLock)
        try
            append!(G,newNodes)
        finally
            unlock(gridLock)
        end
    end
end

function gridHan!(
    G::Grid{T, dim},
    depth::UInt;
    lower::T=-1*one(T),
    upper::T=one(T),
    maxDepth::Union{UInt,Nothing}=nothing,
    offset::Union{Vector{T},Nothing}=nothing,
    nodeCondition=localC,
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
            node = GridNode(G,local_depth,coord;localCondition=nodeCondition)
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
                    local_depth+1;
                    depth=local_depth,
                    scale=(upper-lower)/2
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
    maxDepth::Union{UInt,Nothing}=nothing,
    projectiveCondition=projectiveLocalC
) where {T, dim}
    basisVectors = usualBasis(T, UInt(dim))
    enumeratedCharts = collect(enumerate(PG.charts))
    Threads.@threads for (i,G) in enumeratedCharts
        basisIndex = UInt(((i-1)%dim)+1)
        vectorSign = (-1)^(floor(Int,(i-1)/dim))
        gridHan!(
            G,
            depth;
            lower=lower,
            upper=upper,
            maxDepth=maxDepth,
            offset=basisVectors[basisIndex]*vectorSign,
            nodeCondition=projectiveCondition,
            _split=_projectiveSplitCoordinate(T, basisIndex)
        )
    end
    
    PG.est_condition = maximum(est_condition, PG.charts)
end

function _projectiveSplitCoordinate(::Type{T}, dimExclude::UInt) where T <: Number
    function __projectiveSplitCoordinate(
        v::Vector{T},
        target_depth::UInt;
        depth::UInt = UInt(1),
        scale::T = one(T)
    )
        vProjection = cat(v[1:dimExclude-1],v[dimExclude+1:length(v)];dims=1)
        coords = splitCoordinate(
            vProjection,
            target_depth;
            depth=depth,
            scale=scale
        )
        return [
            cat(u[1:dimExclude-1],v[dimExclude],u[dimExclude:length(u)];dims=1)
            for u in coords
        ]
    end

    return __projectiveSplitCoordinate
end
end #module
