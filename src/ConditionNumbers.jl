module ConditionNumbers
using LinearAlgebra
using ..NormsPolynomials

export localC
export projectiveLocalC

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

    return 1/maximum((scale1, scale2))
end

function projectiveLocalC(
    f::PolynomialSystem{T},
    x::Vector{T};
    image=nothing,
    jacobian=nothing
) where T
    image = isnothing(image) ? f(x) : image
    imageNorm = norm(
        map(((d,y),) -> y/sqrt(d), zip(f.degrees, image)),
        2
    )

    jacobian = isnothing(jacobian) ? f.jacobian(x) : deepcopy(jacobian)
    rowsJacobian = eachrow(jacobian)
    broadcast!(
        (row, d) -> row/d,
        rowsJacobian,
        rowsJacobian,
        f.degrees
    )
    sigmaq = last(svdvals(jacobian))

    return sqrt(length(image)/(imageNorm^2 + sigmaq^2))
end
end # module