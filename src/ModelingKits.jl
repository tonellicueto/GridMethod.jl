## Docs for functions giving compatibility with symbolic modeling packages as HC.ModelKit.jl or ModelingToolkit.jl

"""
    terms(poly; vars, expanded = false)

Returns the list of terms, tuples `(exponent, coefficient)`, of a polynomial `poly`.
The `poly` will be feed by [`get_polys`](@ref).

# Keyword arguments:
 - `vars`: The variables of the polynomial.
 - `expanded`: A boolean indicating whether the polynomial should be expanded before extracting the exponents and coefficients.
"""
function terms end

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
