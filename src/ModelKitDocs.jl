## Docs for functions giving compatibility with symbolic modeling packages as HC.ModelKit.jl or ModelingToolkit.jl

"""
    terms(poly; vars, expanded = false)

Returns the list of terms, tuples `(exponent, coefficient)`, of a polynomial `poly`.
[`get_polys`](@ref) has to return polynomials of the same type as `poly`.

# Keyword arguments:
 - `vars`: The variables of the polynomial.
 - `expanded`: A boolean indicating whether the polynomial should be expanded before extracting the exponents and coefficients.
"""
function terms end

## Getting terms of a system; outputs a vector of iters.
function terms(sys;
               vars = get_vars(sys),
               expanded = true,
               kwargs...)
    return Iterators.map(poly -> terms(poly; vars = vars, expanded = expanded),
                         get_polys(sys))
end


"""
    get_vars(F)

Returns the variables of a system of polynomials `F`.
"""
function get_vars end

"""
    get_polys(F)

Returns the collection `v` of polynomials of a system `F`.
[`terms`](@ref) has to accept polynomials of the element type of `v`.
"""
function get_polys end

"""
    degrees(F)

Returns the degrees of a system of polynomials `F`.
"""
function degrees end
