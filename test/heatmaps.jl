import HomotopyContinuation.ModelKit
const HCMK = ModelKit

using Plots
using GridMethod.GridModule
using GridMethod.GridPlot
using GridMethod.Polynomial
using GridMethod.Han

function main()
    HCMK.@var x,y
    polysys1 = HCMK.System(
        [
        x,
        y,
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
    gridHan!(grid1,UInt(1);maxDepth=UInt(15))
    gridHeatMap(grid1; condition_transform=log)
    savefig("grid1.pdf")

    polysys2 = HCMK.System(
        [
        x - y,
        y,
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

    gridHan!(grid2,UInt(1);maxDepth=UInt(15))
    gridHeatMap(grid2; condition_transform=log)
    savefig("grid2.pdf")

    polysys3 = HCMK.System(
        [
        1000*x - y,
        y,
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

    gridHan!(grid3,UInt(1);maxDepth=UInt(10))
    gridHeatMap(grid3;condition_transform=log)
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

    gridHan!(grid4,UInt(1);maxDepth=UInt(10))
    gridHeatMap(grid4;condition_transform=log)
    savefig("grid4.pdf")

    polysys5 = HCMK.System(
        [
        x^2 + y^2 - 0.5,
        x-0.25,
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

    gridHan!(grid5,UInt(1);maxDepth=UInt(10))
    gridHeatMap(grid5;condition_transform=log)
    savefig("grid5.pdf")

    polysys6 = HCMK.System(
        [
        x^2 + y^2 - 0.5,
        x-(1/sqrt(2.0)),
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
    gridHeatMap(grid6;condition_transform=log)
    savefig("grid6.pdf")
end

main()