module GridPlot
using .Iterators
using ..GridModule
using ..Coordinates
using Plots
using LinearAlgebra

export gridHeatMap

function gridHeatMap(
    grid::Grid{Float64, 2};
    condition_transform=nothing
)
    max_depth = maximum(depth, gridnodes(grid))
    minNode = argmin(coordinates, gridnodes(grid))
    maxNode = argmax(coordinates, gridnodes(grid))

    # First calculate the range of the 
    minCoord = first(sort(splitCoordinate(
        coordinates(minNode),
        max_depth;
        depth=depth(minNode)
    )))
    maxCoord = last(sort(splitCoordinate(
        coordinates(maxNode),
        max_depth;
        depth=depth(maxNode)
    )))

    # Build the x and y ranges for the heatmap
    stepSize = 1/(2^(max_depth-1))
    xRange = collect(range(minCoord[1],maxCoord[1];step=stepSize))
    yRange = collect(range(minCoord[2],maxCoord[2];step=stepSize))

    # Get the "heat" of each coordinate by extracting and sorting
    # the condition numbers of each node.
    points = sort(collect(flatmap(node -> (
            [x[1],x[2],condition(node)]
            for x in splitCoordinate(coordinates(node), max_depth; depth=depth(node))
        ),
        gridnodes(grid)
    )))

    zMatrix = transpose(reshape(
        collect(map(p -> p[3], points)),
        length(xRange),
        length(yRange)
    ))
    if !isnothing(condition_transform)
        broadcast!(condition_transform,zMatrix,zMatrix)
    end

    # Finally, plot
    heatmap(
        xRange,
        yRange,
        zMatrix,
    )
end

end # module