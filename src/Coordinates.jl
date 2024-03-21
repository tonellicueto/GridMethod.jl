module Coordinates
using .Iterators

export splitCoordinate

function findMaxCoordinate(L::Vector{Vector{T}}; dim::UInt=UInt(1)) where T <: Number
    maxV = maxCoordinate(L; dim=dim)
    return findfirst(v -> v == maxV, L)
end

# scale is essentially the side length of the cube centered
# at each coordinate. Hence, to split a cube B centered at
# coordinate x we have to divide scale by 4 to get the center
# of each new cube of side length scale/2.
function splitCoordinate(
    v::Vector{T},
    target_depth::UInt;
    depth::UInt = UInt(1),
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