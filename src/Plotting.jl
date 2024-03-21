
using Plots
pyplot()

using HomotopyContinuation
const HC = HomotopyContinuation
X = @var x,y
X = collect(X)

# using ModelingToolkit
# const MTK = ModelingToolkit
# X = MTK.@variables x,y

# using Revise
using GridMethod
const GM = GridMethod

sys = [GM.rand_poly(X, i) for i in 5:6];
F = HC.System(sys, X);
KF(x...) = GM.K(x, F; expanded = true);
KF(1,2) # first evaluation (requires compilation)

## Adjust parameters to get a nice grid
CC = 28.6
G = GM.grid_Han(KF, CC, 1e-9, 2; depth₀=0);
scattergrid(G)
