## This file contains functions for interact with the ModelingToolkit.jl (MTK) package.

_degree(M::MTK.Symbolic, v::AbstractVector) = [MTK.degree(M, x) for x in v]
function terms(poly::MTK.Num;
                      vars = MTK.get_variables(poly),
                      expanded = false,
                      )
    p = expanded ? MTK.value(poly) : MTK.value(MTK.expand(poly))
    return [(_degree(e, vars), c) for (e, c) in p.dict]
end

get_vars(sys::MTK.NonlinearSystem) = MTK.states(sys)
get_polys(sys::MTK.NonlinearSystem) = [MTK.Num(eq.rhs) for eq in MTK.equations(sys)]
degrees(sys::MTK.NonlinearSystem) = MTK.degree.(get_polys(sys))

function eval_and_J(x, sys::MTK.NonlinearSystem; kwargs...)
    f = MTK.NonlinearFunction(sys;# dvs = states(sys), ps = parameters(sys), u0 = nothing;
                              version = nothing,
                              jac = true,
                              eval_expression = true,
                              sparse = false,
                              simplify = false,
                              kwargs...)
    (f.f(x,[]), f.jac(x,[]))
end
