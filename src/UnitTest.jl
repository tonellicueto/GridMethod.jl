using Test
include("Grid.jl")

# Test suite for GridNode
@testset "GridNode Tests" failfast=true begin
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
@testset "Grid Functionality" begin
    # Create a GridNode instance for testing
    node1 = GridNode{Float64, 3}(2, [1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [1.0 0.0; 0.0 1.0; 0.0 0.0], 10.0)
    node2 = GridNode{Float64, 3}(3, [1.5, 2.5, 3.5], [4.5, 5.5, 6.5], [1.5 0.5; 0.5 1.5; 0.5 0.5], 15.0)

    # Create a Grid instance
    polysys_function = sin # Example polynomial system function
    grid = Grid{Float64, 3, typeof(polysys_function)}(polysys_function, [node1], 5.0)
    
    # Test Grid properties access
    @test polysys(grid) === polysys_function
    @test est_condition(grid) === 5.0
    @test dim(grid) === 3
    @test length(grid) == 1
    @test size(grid) == (1,)
    @test !isempty(grid)
    
    # Test Grid modification
    push!(grid, node2)
    @test length(grid) == 2
    @test grid[2] == node2

    popped_node = pop!(grid)
    @test popped_node == node2
    @test length(grid) == 1

    # Add some more nodes to the Grid
    for i in 0:10
        # Each node is essentially junk
        push!(grid, GridNode{Float64, 3}(
            2,
            [i*1.0, 2.0^i, 3.0^(-1*i)],
            [4.0, 5.0-i, cos(i*6.0)],
            [i*1.0 i^2*1.0; 0.0 1.0; 0.0 0.0],
            1.25*i
        ))
    end 

    # Test iteration over Grid
    for (i, node) in enumerate(grid)
        @test node == grid[i]
    end

    # Test findmax for Grid
    @test condition(findmax(grid)[1]) == 12.5

    @test node2 == GridNode(
        grid,
        3,
        [1.5, 2.5, 3.5],
        [4.5, 5.5, 6.5],
        [1.5 0.5; 0.5 1.5; 0.5 0.5],
        15.0
    )
end