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
    jacobian=nothing,
    jacobian_pinv=nothing
) where T
    image = isnothing(image) ? f(x) : image
    scale1 = norm(
        map(((d,y),) -> y/d, zip(f.degrees, image)),
        Inf
    )

    jacobian = isnothing(jacobian) ? f.jacobian(x) : jacobian
    jacobian_pinv = isnothing(jacobian_pinv) ? pinv(jacobian) : deepcopy(jacobian_pinv)
    cols_pinv = eachcol(jacobian_pinv)
    broadcast!(
        (col, d) -> (d^2)*col,
        cols_pinv,
        cols_pinv,
        f.degrees,
    )
    scale2 = 1/opnorm(jacobian_pinv, Inf)

    return polyNorm1(f)/maximum((scale1, scale2))
end
end # module