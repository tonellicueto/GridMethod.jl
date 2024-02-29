import HomotopyContinuation.ModelKit
const HCMK = ModelKit

using .Iterators
using Test
using LinearAlgebra
using GridMethod.Polynomial
using GridMethod.Norms
using GridMethod.ConditionNumbers
const CN = ConditionNumbers
using GridMethod.GridModule
using GridMethod.GridBatch
using GridMethod.Han


# Test suite for GridNode
@testset "GridNode Basic Tests" failfast=true begin
    # Create an instance of GridNode for testing
    coordinates_ = [1.0, 2.0, 3.0]
    image_ = [4.0, 5.0, 6.0]
    jacobian_ = [1.0 0.0; 0.0 1.0; 0.0 0.0]
    condition_ = 10.0
    node = GridNode{Float64, 3}(2, coordinates_, image_, jacobian_, condition_)

    # Test the depth function
    @test depth(node) == 2

    # Test the coordinates function
    @test coordinates(node) == [1.0, 2.0, 3.0]

    # Test the image function
    @test image(node) == [4.0, 5.0, 6.0]

    # Test the jacobian function
    @test jacobian(node) == [1.0 0.0; 0.0 1.0; 0.0 0.0]

    # Test the condition function
    @test condition(node) == 10.0

    # Test the dim function
    @test dim(node) == 3

    coordinates2 = [2.0, 3.0, 4.0]
    image2 = [5.0, 6.0, 7.0]
    jacobian2 = [2.0 -10.0; 5.0 2.0; 1.0 0.0]
    condition2 = 1000.0
    node2 = GridNode{Float64, 3}(2, coordinates2, image2, jacobian2, condition2)

    @test node < node2
    @test node == node
    @test node != node2

end

# Start of the test suite for Grid
@testset "Grid Basic Functionality" failfast=true begin
    # Create a Grid instance
    HCMK.@var x,y,z
    polysys_ = HCMK.System(
        [x-y, y-z, x^2 -z^3];
        variables=[x,y,z]
    )
    jacobian_ = v -> HCMK.jacobian(polysys_, v)

    gridPolySys::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys_(v),
        jacobian_,
        HCMK.degrees(polysys_),
        HCMK.support_coefficients(polysys_)[2]
    )
    # Make a Grid for testing
    grid = Grid{Float64, 3}(gridPolySys, [], nothing)
    
    # Test Grid properties access
    @test polysys(grid) === gridPolySys
    @test est_condition(grid) === nothing 
    @test dim(grid) == 3
    @test length(grid) == 0
    @test size(grid) == (0,)
    @test isempty(grid)

    # Create GridNode instances for testing
    coordinate1 = [1.0, 2.0, 3.0]
    node1 = GridNode(
        grid,
        2,
        coordinate1
    )
    coordinate2 = [1.5, 2.5, 3.5]
    node2 = GridNode(
        grid,
        3,
        coordinate2,
    )

    # Test Grid modification
    push!(grid, node1)
    push!(grid, node2)
    @test length(grid) == 2
    @test grid[2] == node2

    popped_node = pop!(grid)
    @test popped_node == node2
    @test length(grid) == 1
    push!(grid, node2)

    # Test iteration over Grid
    for (i, node) in enumerate(grid)
        @test node == grid[i]
    end

    # Test node evaluation
    for _ in 1:length(grid)
        node = pop!(grid)
        newNode = GridNodeEvaluate(grid, node)
        pushfirst!(grid, newNode)
    end

    for node in grid
        @test node.image == polysys(grid)(coordinates(node))
        @test node.jacobian == polysys(grid).jacobian(coordinates(node))
        @test node.condition == localC(polysys(grid), coordinates(node))
    end
end

@testset "Norm Tests" failfast=true begin
    HCMK.@var x,y,z
    polysys_ = HCMK.System(
        [x-y, x+y-z, x^2 -z^3];
        variables=[x,y,z]
    )
    jacobian_ = v -> HCMK.jacobian(polysys_, v)

    polysystem::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys_(v),
        jacobian_,
        HCMK.degrees(polysys_),
        HCMK.support_coefficients(polysys_)[2]
    )
    
    @test polyNorm1(polysystem) == 3.0 

    # hinfnorm
    @test hinfNorm([0.1, -0.2]) == 1.0
    @test hinfNorm([0.1, -2]) == 2.0

    # matrixInfPNorm
    A = [1 0 1; 1 1 0]
    @test matrixInfPNorm(A; p=2)==sqrt(8)
end

@testset "Condition Number Tests" failfast=true begin
    @test CN._vector_power(2.0, [0.0, -1.0, 1.5]) == [1.0, 0.5, 2*sqrt(2)]

    # Set up polynomial system for condition number testing
    HCMK.@var x1,x2,x3
    polysys_ = HCMK.System(
        [x1+x2-x3, x1^2 -x3^3];
        variables=[x1,x2,x3]
    )
    jacobian_ = v -> HCMK.jacobian(polysys_, v)

    polysystem::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys_(v),
        jacobian_,
        HCMK.degrees(polysys_),
        HCMK.support_coefficients(polysys_)[2]
    )

    testVectors = [
        [1.0,1.0,1.0],
        [0.0,0.25,-0.5],
        [0.375,-0.1,-0.8],
        [0.0,0.0,0.0],
    ]
    for v in testVectors
        # Work out condition number manually
        d = polysystem.degrees
        Δ = Diagonal(d)
        vhinf = hinfNorm(v)

        A1 = inv(Δ)*Diagonal(CN._vector_power(vhinf, -1.0*d))
        scale1 = norm(A1*polysystem(v), Inf)

        jacobianV = polysystem.jacobian(v)
        A2 = (
            pinv(jacobianV)
            *Diagonal(CN._vector_power(vhinf, d-ones(Float64, size(d)[1])))
            *Δ^2
        )
        scale2 = 1/opnorm(A2, Inf)

        condV = polyNorm1(polysystem)/maximum([scale1,scale2])

        @test localC(polysystem, v) == condV
    end
end

# Test coordinate batching
@testset "GridBatch " failfast=true begin
    testBatches = sort(collect(flatten(
        GridBatch.batchGrid(Float64,3,2,4,-1.0,1.0)
    )))

    expectedBatches = sort([
        [l,r]
        for l in range(-0.75,0.75,step=0.5)
        for r in range(-0.75,0.75,step=0.5)
    ])

    @test testBatches == expectedBatches
end

# Test coordinate batching
@testset "Han test" failfast=true begin
    # Create a Grid instance
    HCMK.@var x,y,z
    polysys_ = HCMK.System(
        [x-y, y-z, x^2 -z^3];
        variables=[x,y,z]
    )
    jacobian_ = v -> HCMK.jacobian(polysys_, v)

    gridPolySys::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys_(v),
        jacobian_,
        HCMK.degrees(polysys_),
        HCMK.support_coefficients(polysys_)[2]
    )
    # Make a Grid for testing
    grid = Grid{Float64, 3}(gridPolySys, [], nothing)

    @test length(grid) == 0
    gridHan!(grid,5)
    @test length(grid) > 0
end