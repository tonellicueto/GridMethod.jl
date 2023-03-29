
using PlotlyJS

using HomotopyContinuation
const HC = HomotopyContinuation

X = @var x,y
X = collect(X)

# using ModelingToolkit
# const MTK = ModelingToolkit

# using Revise
using GridMethod
const GM = GridMethod

sys = [GM.rand_poly(X, i) for i in 5:6];
F = HC.System(sys, X);
KF(x,y) = GM.K([x,y], F; expanded = true);
KF(1,2) # first evaluation (requires computation)

## Adjust parameters to get a nice grid
CC = 131
G = GM.refine_grid(KF, CC, 1e-9, 2; depth₀=0);
P = GM.points(G, CC)

points = [p[2] for p in P];
x = first.(points);
y = last.(points);


trace1 = histogram2dcontour(
    x = x,
    y = y,
    colorscale = "Blues",
    reversescale = true,
    xaxis = "x",
    yaxis = "y"
)

trace2 = scatter(
    x = x,
    y = y,
    xaxis = "x",
    yaxis = "y",
    mode = "markers",
    marker = attr(
        color = "rgba(0,0,0,0.7)",
        size = 5
    )
)

trace3 = histogram(
    y = y,
    xaxis = "x2",
    marker = attr(
        color = "rgba(0,0,0,1)"
    )
)

trace4 = histogram(
    x = x,
    yaxis = "y2",
    marker = attr(
        color = "rgba(0,0,0,1)"
    )
)


layout = Layout(
    autosize = false,
    xaxis = attr(
        zeroline = false,
        domain = [0,0.85],
        showgrid = false
    ),

    yaxis = attr(
        zeroline = false,
        domain = [0,0.85],
        showgrid = false
    ),
    xaxis2 = attr(
        zeroline = false,
        domain = [0.85,1],
        showgrid = false
    ),
    yaxis2 = attr(
        zeroline = false,
        domain = [0.85,1],
        showgrid = false
    ),
    height = 600,
    width = 600,
    bargap = 0,
    hovermode = "closest",
    showlegend = false
)

plot([trace1, trace2, trace3, trace4], layout)
