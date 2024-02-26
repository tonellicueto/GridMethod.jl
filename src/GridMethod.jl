module GridMethod

#Indicate modules for use
using UnPack # @unpack, @pack!
using DocStringExtensions: SIGNATURES, TYPEDEF

using Combinatorics
import LinearAlgebra as LA

using DocStringExtensions

#Include files with extracode
include("Polynomial.jl")
using .Polynomial

#export rand_poly
#include("Utils.jl")

#export Onorm, Ocond, Wnorm, Wcond
include("Norms.jl")
using .Norms

include("ConditionNumbers.jl")
using .ConditionNumbers

#using RecipesBase # For Plots recipes
#include("Plots.jl")

# using StaticArrays # Maybe for GridNode
include("Grid.jl")
using .GridModule

#include("Han.jl")
#
#using Requires
#
#include("ModelKitDocs.jl")
#function __init__()
#    @require HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327" begin
#        import .HomotopyContinuation.ModelKit as MK
#        include("HC.ModelKit.jl")
#    end
#    @require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" begin
#        import .ModelingToolkit as MTK
#        include("ModelingToolkit.jl")
#    end
#end


end
