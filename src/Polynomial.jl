module Polynomial
export PolynomialSystem

struct PolynomialSystem{T}
    evaluate::Function
    jacobian::Function
    degrees::Vector{Int}
    coefficients::Vector{Vector{T}}
end

function (p::PolynomialSystem{T})(x::Vector{T}) where T
    p.evaluate(x)
end
end