module Coordinates
using .Iterators

export splitCoordinate

function splitCoordinate(
    v::Vector{T},
    target_depth::UInt;
    depth::UInt = UInt(0),
    scale::T = one(T)
) where {T}
    if depth == target_depth
        return [v]
    end
    new_depth = depth + 1
    vSplit = [
        [w_/2^new_depth for w_ in w] + v
        for w in product(repeated([-one(T)*scale,one(T)*scale],length(v))...)
    ]
    return collect(flatmap(
        u -> splitCoordinate(u,target_depth;depth=new_depth,scale=scale),
        vSplit
    ))
end
end # module