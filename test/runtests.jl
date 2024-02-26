import HomotopyContinuation.ModelKit
const HCMK = ModelKit

using Test
using GridMethod.Polynomial
using GridMethod.Norms
using GridMethod.ConditionNumbers
using GridMethod.GridModule


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

    # TODO Test findmax for Grid
    # @test condition(findmax(grid)[1]) == 12.5
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
    @test ConditionNumbers._vector_power(2.0, [0.0, -1.0, 1.5]) == [1.0, 0.5, 2*sqrt(2)]

    # Set up polynomial system for condition number testing
    HCMK.@var x,y,z
    polysys_ = HCMK.System(
        [x+y-z, x^2 -z^3];
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

    _testVal = localC(polysystem, [1.0,1.0,1.0])
    @info "$_testVal"
    @test !isnothing(_testVal)
end