module GridBatch
using .Iterators
export batchGrid

function batchGrid(
    ::Type{T},
    depth::Integer,
    dim::Integer,
    numBatches::Integer,
    lower::T,
    upper::T
) where T <: Number
    scale = (upper - lower)/2^depth
    offset::Vector{T} = repeat([lower*one(T)], dim)

    batchSize = (2^(depth-1))÷numBatches
    batches = [(i,batchSize) for i in 0:numBatches-1]
    if (2^(depth-1))%numBatches != 0
        push!(batches,(numBatches+1,(2^(depth-1))%numBatches))
    end

    return [
        _makeBatch(T, b...,offset,dim,depth,scale)
        for b in batches
    ]
end

function _makeBatch(
    ::Type{T},
    batchNumber::Integer,
    batchSize::Integer,
    offset::Vector{T},
    dim::Integer,
    depth::Integer,
    scale::T
) where T <: Number
    gridRange = range(1,2^depth-1,step=2)
    rangeStart = batchNumber*batchSize
    batchRange = map(
        i -> 2*i + 1,
        range(rangeStart,rangeStart+batchSize-1)
    )

    if dim > 1
        gridTuples = product(
            batchRange,
            repeated(gridRange, dim-1)...
        )
    else
        gridTuples = ((i,) for i in batchRange)
    end

    return map(
        t -> offset + [n*scale for n in t],
        gridTuples
    )
end
end #module