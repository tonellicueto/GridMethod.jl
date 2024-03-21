import HomotopyContinuation.ModelKit
const HCMK = ModelKit

using Plots
using GridMethod.GridModule
using GridMethod.GridPlot
using GridMethod.Polynomial
using GridMethod.Han

function grid_groupbydepth(G::Grid{T, dim}) where {T, dim}
    I = unique(sort(depth.(G)))
    out = Dict(I .=> [empty(G) for _ in I ])
    for g in G
        d = depth(g)
        push!(out[d], g)
    end
    return out
end

set_forplot(G) = eachrow(reduce(hcat, map(coordinates, gridnodes(G))))

@recipe function f(G::Grid{Float64,2};
                   r = 2,
                   R = 20,
                   f_ms = i -> r + R*(.5)^i, # Marker size
                   f_ma = (d, dmin, dmax) -> (d - dmin)/(dmax - dmin), # Mark alpha/opacity
                   f_msw = (d, dmin, dmax) -> r/10 + (1 - f_ma(d, dmin, dmax))/r, # Border width (0 = No border)
                   color = "rgb(238,37,35)"
                   )
    Dict_depth = grid_groupbydepth(G)
    I = sort!(collect(keys(Dict_depth))) # Sort it for the legend
    dmin = first(I)
    dmax = last(I)
    # Set default values
    axis --> ([-1.1,1.1],)
    # xaxis --> ("x",)
    # yaxis --> ("y",)
    legend --> :outerright
    legendfontsize --> 10
    thickness_scaling --> 1
    shift = dmin -1
    # Plot each depth group
    for d in I
        x, y = set_forplot(Dict_depth[d])
        @series begin
            seriestype := :scatter
            label := "Depth $d"
            # markershape := :circle
            ms := f_ms(d - shift) # Mark size
            mc := color # Mark color
            ma := f_ma(d, dmin, dmax) # Mark alpha/opacity
            msw := f_msw(d, dmin, dmax) # Border/stroke width
            # msw := 0 # No border/stroke
            msc := color # Border color
            msa := f_ma(d, dmin, dmax)/r # Border alpha
            x, y
        end
    end
end

function main()
    HCMK.@var x,y
    polysys1 = HCMK.System(
        [
        y,
        x
        ];
        variables=[x,y]
    )
    jacobian1 = v -> HCMK.jacobian(polysys1, v)

    gridPolySys1::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys1(v),
        jacobian1,
        HCMK.degrees(polysys1),
        HCMK.support_coefficients(polysys1)[2]
    )
    # Make a Grid for testing
    grid1 = Grid{Float64, 2}(gridPolySys1, [], nothing)
    gridHan!(grid1,UInt(1);maxDepth=UInt(12))
    plot(grid1)
    savefig("grid1.pdf")

    polysys2 = HCMK.System(
        [
        y - x,
        x
        ];
        variables=[x,y]
    )
    jacobian2 = v -> HCMK.jacobian(polysys2, v)

    gridPolySys2::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys2(v),
        jacobian2,
        HCMK.degrees(polysys2),
        HCMK.support_coefficients(polysys2)[2]
    )
    # Make a Grid for testing
    grid2 = Grid{Float64, 2}(gridPolySys2, [], nothing)

    gridHan!(grid2,UInt(1);maxDepth=UInt(12))
    plot(grid2)
    savefig("grid2.pdf")

    polysys3 = HCMK.System(
        [
        1000*y - x,
        x
        ];
        variables=[x,y]
    )
    jacobian3 = v -> HCMK.jacobian(polysys3, v)

    gridPolySys3::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys3(v),
        jacobian3,
        HCMK.degrees(polysys3),
        HCMK.support_coefficients(polysys3)[2]
    )
    # Make a Grid for testing
    grid3 = Grid{Float64, 2}(gridPolySys3, [], nothing)

    gridHan!(grid3,UInt(1);maxDepth=UInt(12))
    plot(grid3)
    savefig("grid3.pdf")

    polysys4 = HCMK.System(
        [
        x^2 + y^2 - 0.5,
        x,
        ];
        variables=[x,y]
    )
    jacobian4 = v -> HCMK.jacobian(polysys4, v)

    gridPolySys4::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys4(v),
        jacobian4,
        HCMK.degrees(polysys4),
        HCMK.support_coefficients(polysys4)[2]
    )
    # Make a Grid for testing
    grid4 = Grid{Float64, 2}(gridPolySys4, [], nothing)

    gridHan!(grid4,UInt(1);maxDepth=UInt(12))
    plot(grid4)
    savefig("grid4.pdf")

    polysys5 = HCMK.System(
        [
        x^2 + y^2 - 0.5,
        sqrt(2.0)*x-1.0,
        ];
        variables=[x,y]
    )
    jacobian5 = v -> HCMK.jacobian(polysys5, v)

    gridPolySys5::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys5(v),
        jacobian5,
        HCMK.degrees(polysys5),
        HCMK.support_coefficients(polysys5)[2]
    )
    # Make a Grid for testing
    grid5 = Grid{Float64, 2}(gridPolySys5, [], nothing)

    gridHan!(grid5,UInt(1);maxDepth=UInt(12))
    plot(grid5)
    savefig("grid5.pdf")

    polysys6 = HCMK.System(
        [
        x^2 + y^2 - 0.5,
        sqrt(2.0)*x-1.0,
        ];
        variables=[x,y]
    )
    jacobian6 = v -> HCMK.jacobian(polysys6, v)

    gridPolySys6::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys6(v),
        jacobian6,
        HCMK.degrees(polysys6),
        HCMK.support_coefficients(polysys6)[2]
    )
    # Make a Grid for testing
    grid6 = Grid{Float64, 2}(gridPolySys6, [], nothing)

    gridHan!(grid6,UInt(1);maxDepth=UInt(10))
    plot(grid6)
    savefig("grid6.pdf")
end

main()