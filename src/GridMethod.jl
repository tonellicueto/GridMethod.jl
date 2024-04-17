module GridMethod

include("NormsPolynomials.jl")
using .NormsPolynomials

include("ConditionNumbers.jl")
using .ConditionNumbers

include("Grid.jl")
using .GridModule

include("Coordinates.jl")
using .Coordinates

include("Han.jl")
using .Han
end
