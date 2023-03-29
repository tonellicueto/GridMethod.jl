module GridMethod

#Indicate modules for use
using UnPack # @unpack, @pack!
# using DocStringExtensions: SIGNATURES, TYPEDEF
# using RecipesBase # For Plots recipes
# using InvertedIndices

using Combinatorics
import LinearAlgebra as LA

using DocStringExtensions

#Include files with extracode
export rand_poly
include("Utils.jl")

export norm1, C, Wnorm, K
include("Norms.jl")
include("Const.jl")

using StaticArrays
include("Grid_Refine.jl")
include("Han.jl")

using Requires

function __init__()
    @require HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327" begin
        # import .HomotopyContinuation as HC
        import .HomotopyContinuation.ModelKit as MK
        include("ModelingKits.jl")
        include("HC.ModelKit.jl")
    end
    @require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" begin
        @eval import Symbolics
        import .ModelingToolkit as MTK
        include("ModelingKits.jl")
        include("ModelingToolkit.jl")
    end
end


end
