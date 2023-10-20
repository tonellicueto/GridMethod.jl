
import HomotopyContinuation as HC

X = @var x,y
X = collect(X)
degs = [7, 9]

# using Revise
using GridMethod
sys = [GridMethod.rand_poly(X, d; homo = false) for d in degs];

F = HC.System(sys, X);

OcondF(x,y) = inv(GridMethod.Ocond(F, [x,y]; compiled = true));
OcondF(.1,.2) # Compile OcondF

C = maximum(degs)*2.;
# HanF = GridMethod.grid_Han(OcondF, C, 1e-9, 2; depth₀=0);
HanF = GridMethod.grid_Han(OcondF, C, 1e-9, 2);

using Plots
# Do not execute:
# ENV["QT_QPA_PLATFORM"] = :wayland
# pyplot() # Use Python as backend

scattergrid(HanF)
