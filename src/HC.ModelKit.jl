## This file contains functions for interact with the HomotopyContinuation.ModelKit.jl (MK) package.

function terms(poly::MK.Expression;
               vars = MK.variables(poly),
               expanded = false,
               )
    p = expanded ? poly : MK.expand(poly)
    return [(e, MK.to_number(c)) for (e, c) in MK.to_dict(p, vars)]
end

get_vars(F::MK.System) = MK.variables(F)
get_polys(F::MK.System) = MK.expressions(F)
degrees(F::MK.System) = MK.degrees(F)

function eval_and_J(x, F::MK.System; compiled = true, kwargs...)
    f = compiled ? MK.CompiledSystem(F; optimizations=true) : MK.InterpretedSystem(F; optimizations=true)
    u = Vector{Any}(undef, size(F, 1))
    U = Matrix{Any}(undef, size(F))
    MK.evaluate_and_jacobian!(u, U, f, x, nothing)
    fx = MK.to_smallest_eltype(u)
    Jx = MK.to_smallest_eltype(U)
    (x, fx, Jx)
end
