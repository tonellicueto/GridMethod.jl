module GridMethod

#Include files with extracode
include("Polynomial.jl")
using .Polynomial

include("Norms.jl")
using .Norms

include("ConditionNumbers.jl")
using .ConditionNumbers

include("Grid.jl")
using .GridModule

include("Coordinates.jl")
using .Coordinates

include("Han.jl")
using .Han
end
