module ConditionNumbers
using LinearAlgebra
using ..Polynomial
using ..Norms

export localC

function _vector_power(x::T, v::Vector{T}) where T <: Number
    return map(y -> x^y, v)
end

function localC(
    f::PolynomialSystem{T},
    x::Vector{T};
    image=nothing,
    jacobian=nothing
) where T
    image = isnothing(image) ? f(x) : image
    scale1 = norm(
        map(((d,y),) -> y/d, zip(f.degrees, image)),
        Inf
    )

    jacobian = isnothing(jacobian) ? f.jacobian(x) : deepcopy(jacobian)
    rowsJacobian = eachrow(jacobian)
    broadcast!(
        (row, d) -> row/d^2,
        rowsJacobian,
        rowsJacobian,
        f.degrees
    )
    scale2 = last(filter(x -> x!=0.0, svdvals(jacobian)))

    return polyNorm1(f)/maximum((scale1, scale2))
end
end # module