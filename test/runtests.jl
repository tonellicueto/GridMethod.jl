import HomotopyContinuation.ModelKit
const HCMK = ModelKit

using .Iterators
using Test
using LinearAlgebra
using GridMethod.NormsPolynomials
using GridMethod.ConditionNumbers
const CN = ConditionNumbers
using GridMethod.GridModule
using GridMethod.Han
using GridMethod.Coordinates

function _convertHCMKDegreesToPolynomial(degrees::Vector{Matrix{T}}) where T
    return [
        collect(eachcol(M))
        for M in degrees
    ]
end

@testset "Convert support coefficients to Polynomial format" failfast=true begin
    HCMK.@var x,y,z
    P = HCMK.System(
        [x-y, y-z, x^2 -z^3, x^2*y^3-2.0*z^5+x*y^2*z^2];
        variables=[x,y,z]
    )

    L = _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(P)[1])
    @test L[1][1] == [1, 0, 0] 
    @test L[1][2] == [0, 1, 0] 
    @test L[2][1] == [0, 1, 0] 
    @test L[2][2] == [0, 0, 1] 
    @test L[3][1] == [0, 0, 3] 
    @test L[3][2] == [2, 0, 0] 
    @test L[4][1] == [2, 3, 0] 
    @test L[4][2] == [1, 2, 2] 
    @test L[4][3] == [0, 0, 5] 
end


# Test suite for GridNode
@testset "GridNode Basic Tests" failfast=true begin
    # Create an instance of GridNode for testing
    coordinates_ = [1.0, 2.0, 3.0]
    image_ = [4.0, 5.0, 6.0]
    jacobian_ = [1.0 0.0; 0.0 1.0; 0.0 0.0]
    condition_ = 10.0
    node = GridNode{Float64, 3}(UInt(2), coordinates_, image_, jacobian_, condition_)

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
    node2 = GridNode{Float64, 3}(UInt(2), coordinates2, image2, jacobian2, condition2)

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
        HCMK.support_coefficients(polysys_)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys_)[1])
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
        UInt(2),
        coordinate1
    )
    coordinate2 = [1.5, 2.5, 3.5]
    node2 = GridNode(
        grid,
        UInt(3),
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
        newNode = GridNode(grid, depth(node), coordinates(node))
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
        HCMK.support_coefficients(polysys_)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys_)[1])
    )
    
    #@test polyNorm1(polysystem) == 3.0 

    # hinfnorm
    @test hinfNorm([0.1, -0.2]) == 1.0
    @test hinfNorm([0.1, -2]) == 2.0

    # matrixInfPNorm
    A = [1 0 1; 1 1 0]
    @test matrixInfPNorm(A; p=2)==sqrt(8)

    degree = 3
    monomials = [[3, 0, 0], [1, 1, 1,], [1,0,0], [0, 1, 2]]
    coefficients = [-1.0, 2.0, -3.0, 10.0]
    @test weylNorm(degree, monomials, coefficients)==sqrt(379.0)
end

@testset "Condition Number Tests" failfast=true begin
    @test CN._vector_power(2.0, [0.0, -1.0, 1.5]) == [1.0, 0.5, 2*sqrt(2)]

    # Set up polynomial system for condition number testing
    HCMK.@var x1,x2
    polysys_ = HCMK.System(
        [1000*x1-x2,x2];
        variables=[x1,x2]
    )
    jacobian_ = v -> HCMK.jacobian(polysys_, v)

    polysystem::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys_(v),
        jacobian_,
        HCMK.degrees(polysys_),
        HCMK.support_coefficients(polysys_)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys_)[1])
    )

    testVectors = [
        [1.0,1.0],
        [0.0,0.25],
        [0.375,-0.1],
        [0.0,0.0],
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
            pinv(
                Diagonal(CN._vector_power(vhinf, d-ones(Float64, size(d)[1])))
                *Δ^2
            )
            *jacobianV
        ) 
        scale2 = last(filter(x -> x!=0.0, svdvals(A2))) 

#        TODO Failing
#        condV = PolyNorm1(polysystem)/maximum([scale1,scale2])
#
#        @test localC(polysystem, v)==condV
    end

    # TODO finish unit testing projectiveLocalC
    polysys2_ = HCMK.System(
        [x1^2*x2 - x1*x2^2];
        variables=[x1,x2]
    )
    jacobian2_ = v -> HCMK.jacobian(polysys2_, v)

    polysystem2::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys2_(v),
        jacobian2_,
        HCMK.degrees(polysys2_),
        HCMK.support_coefficients(polysys2_)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys2_)[1])
    )

end

@testset "Han test" failfast=true begin
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
        HCMK.support_coefficients(polysys1)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys1)[1])
    )
    # Make a Grid for testing
    grid1 = Grid{Float64, 2}(gridPolySys1, [], nothing)

    @test length(grid1) == 0
    gridHan!(grid1,UInt(1);maxDepth=UInt(15))
    @test length(grid1) > 0

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
        HCMK.support_coefficients(polysys2)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys2)[1])
    )
    # Make a Grid for testing
    grid2 = Grid{Float64, 2}(gridPolySys2, [], nothing)

    @test length(grid2) == 0
    gridHan!(grid2,UInt(1);maxDepth=UInt(15))
    @test length(grid2) > 0

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
        HCMK.support_coefficients(polysys3)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys3)[1])
    )
    # Make a Grid for testing
    grid3 = Grid{Float64, 2}(gridPolySys3, [], nothing)

    @test length(grid3) == 0
    gridHan!(grid3,UInt(1);maxDepth=UInt(10))
    @test length(grid3) > 0

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
        HCMK.support_coefficients(polysys4)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys4)[1])
    )
    # Make a Grid for testing
    grid4 = Grid{Float64, 2}(gridPolySys4, [], nothing)

    @test length(grid4) == 0
    gridHan!(grid4,UInt(1);maxDepth=UInt(10))
    @test length(grid4) > 0

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
        HCMK.support_coefficients(polysys5)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys5)[1])
    )
    # Make a Grid for testing
    grid5 = Grid{Float64, 2}(gridPolySys5, [], nothing)

    @test length(grid5) == 0
    gridHan!(grid5,UInt(1);maxDepth=UInt(10))
    @test length(grid5) > 0

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
        HCMK.support_coefficients(polysys6)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys6)[1])
    )
    # Make a Grid for testing
    grid6 = Grid{Float64, 2}(gridPolySys6, [], nothing)

    @test length(grid6) == 0
    gridHan!(grid6,UInt(1);maxDepth=UInt(10))
    @test length(grid6) > 0

    # Test min depth increase
    minDepthGrid6 = minimum(depth, gridnodes(grid6))
    maxDepthGrid6 = maximum(depth, gridnodes(grid6))

    increaseMinDepth!(grid6, UInt(5);scale=2.0)
    newMinDepthGrid6 = minimum(depth, gridnodes(grid6))
    newMaxDepthGrid6 = maximum(depth, gridnodes(grid6))
    @test newMinDepthGrid6==max(5,minDepthGrid6)
    @test newMaxDepthGrid6==maxDepthGrid6

    oldMinDepthGrid6 = minimum(depth, gridnodes(grid6))
    oldMaxDepthGrid6 = maximum(depth, gridnodes(grid6))
    increaseDepth!(grid6, UInt(2);scale=2.0)
    newMinDepthGrid6 = minimum(depth, gridnodes(grid6))
    newMaxDepthGrid6 = maximum(depth, gridnodes(grid6))
    @test newMinDepthGrid6==oldMinDepthGrid6+2
    @test newMaxDepthGrid6==oldMaxDepthGrid6+2

    oldMinDepthGrid6 = minimum(depth, gridnodes(grid6))
    oldMaxDepthGrid6 = maximum(depth, gridnodes(grid6))
    increaseMinDepth!(grid6, UInt(8);scale=2.0,nodeFilter=(node)->depth(node)≤10)
    newMinDepthGrid6 = minimum(depth, gridnodes(grid6))
    newMaxDepthGrid6 = maximum(depth, gridnodes(grid6))
    @test newMinDepthGrid6==max(oldMinDepthGrid6,8)
    @test newMaxDepthGrid6≤10

    oldMinDepthGrid6 = minimum(depth, gridnodes(grid6))
    oldMaxDepthGrid6 = maximum(depth, gridnodes(grid6))
    increaseDepth!(grid6, UInt(1);scale=2.0,nodeFilter=(node)->depth(node)≤10)
    newMinDepthGrid6 = minimum(depth, gridnodes(grid6))
    newMaxDepthGrid6 = maximum(depth, gridnodes(grid6))
    @test newMinDepthGrid6==max(oldMinDepthGrid6,9)
    @test newMaxDepthGrid6≤10
end

@testset "Coordinates test" failfast=true begin
    x = [0.0,0.0,0.0]
    splitz = splitCoordinate(x, UInt(1))
    @test length(splitz) == 1
    @test sort(splitz) == [[0.0,0.0,0.0]]

    splitx = splitCoordinate(x, UInt(2);scale=2.0)
    @test length(splitx) == 8
    @test sort(splitx) == [
        [-0.5,-0.5,-0.5],
        [-0.5,-0.5,0.5],
        [-0.5,0.5,-0.5],
        [-0.5,0.5,0.5],
        [0.5,-0.5,-0.5],
        [0.5,-0.5,0.5],
        [0.5,0.5,-0.5],
        [0.5,0.5,0.5],
    ]

    y = [0.0,0.0]
    splity = splitCoordinate(y, UInt(3);scale=0.5)
    @test length(splity) == 16
    @test sort(splity) == [
        [-0.1875,-0.1875],
        [-0.1875,-0.0625],
        [-0.1875,0.0625],
        [-0.1875,0.1875],
        [-0.0625,-0.1875],
        [-0.0625,-0.0625],
        [-0.0625,0.0625],
        [-0.0625,0.1875],
        [0.0625,-0.1875],
        [0.0625,-0.0625],
        [0.0625,0.0625],
        [0.0625,0.1875],
        [0.1875,-0.1875],
        [0.1875,-0.0625],
        [0.1875,0.0625],
        [0.1875,0.1875],
    ]

    basisR2 = usualBasis(Float64, UInt(2))
    @test length(basisR2) == 2
    @test basisR2[1] == [1.0, 0.0]
    @test basisR2[2] == [0.0, 1.0]
end

@testset "ProjectiveGrid tests" failfast=true begin
    HCMK.@var x,y
    polysys1 = HCMK.System(
        [
        x + y,
        x - y,
        ];
        variables=[x,y]
    )
    jacobian1 = v -> HCMK.jacobian(polysys1, v)

    gridPolySys1::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys1(v),
        jacobian1,
        HCMK.degrees(polysys1),
        HCMK.support_coefficients(polysys1)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys1)[1])
    )
    PG1::ProjectiveGrid{Float64, UInt(2)} = ProjectiveGrid(
        gridPolySys1,
        UInt(2)
    )

    @test size(PG1.charts)[1] == 4
    @test isnothing(PG1.est_condition)

    projectiveGridHan!(PG1,UInt(1);maxDepth=UInt(10))

    HCMK.@var z
    polysys2 = HCMK.System(
        [
        x + y + z,
        ];
        variables=[x,y,z]
    )
    jacobian2 = v -> HCMK.jacobian(polysys2, v)

    gridPolySys2::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys2(v),
        jacobian2,
        HCMK.degrees(polysys2),
        HCMK.support_coefficients(polysys2)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys2)[1])
    )
    PG2::ProjectiveGrid{Float64, UInt(3)} = ProjectiveGrid(
        gridPolySys2,
        UInt(3)
    )

    @test size(PG2.charts)[1] == 6
    @test isnothing(PG2.est_condition)

    projectiveGridHan!(PG2,UInt(1);maxDepth=UInt(10))

    polysys3 = HCMK.System(
        [
        x + y + z,
        x + y + (1+1/100)*z,
        ];
        variables=[x,y,z]
    )
    jacobian3 = v -> HCMK.jacobian(polysys3, v)

    gridPolySys3::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys3(v),
        jacobian3,
        HCMK.degrees(polysys3),
        HCMK.support_coefficients(polysys3)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys3)[1])
    )
    PG3::ProjectiveGrid{Float64, UInt(3)} = ProjectiveGrid(
        gridPolySys3,
        UInt(3)
    )

    @test size(PG3.charts)[1] == 6
    @test isnothing(PG3.est_condition)

    projectiveGridHan!(PG3,UInt(1);maxDepth=UInt(10))

    polysys4 = HCMK.System(
        [
        x^2 + y^2 - z^2,
        z,
        ];
        variables=[x,y,z]
    )
    jacobian4 = v -> HCMK.jacobian(polysys4, v)

    gridPolySys4::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys4(v),
        jacobian4,
        HCMK.degrees(polysys4),
        HCMK.support_coefficients(polysys4)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys4)[1])
    )
    PG4::ProjectiveGrid{Float64, UInt(3)} = ProjectiveGrid(
        gridPolySys4,
        UInt(3)
    )

    @test size(PG4.charts)[1] == 6
    @test isnothing(PG4.est_condition)

    projectiveGridHan!(PG4,UInt(1);maxDepth=UInt(10))

    polysys5 = HCMK.System(
        [
        x^2 - y*z,
        z,
        ];
        variables=[x,y,z]
    )
    jacobian5 = v -> HCMK.jacobian(polysys5, v)

    gridPolySys5::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys5(v),
        jacobian5,
        HCMK.degrees(polysys5),
        HCMK.support_coefficients(polysys5)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys5)[1])
    )
    PG5::ProjectiveGrid{Float64, UInt(3)} = ProjectiveGrid(
        gridPolySys5,
        UInt(3)
    )

    @test size(PG5.charts)[1] == 6
    @test isnothing(PG5.est_condition)

    projectiveGridHan!(PG5,UInt(1);maxDepth=UInt(10))

    polysys6 = HCMK.System(
        [
        x^2 - y^2 - 1,
        z,
        ];
        variables=[x,y,z]
    )
    jacobian6 = v -> HCMK.jacobian(polysys6, v)

    gridPolySys6::PolynomialSystem{Float64} =
    PolynomialSystem{Float64}(
        v -> polysys6(v),
        jacobian6,
        HCMK.degrees(polysys6),
        HCMK.support_coefficients(polysys6)[2],
        _convertHCMKDegreesToPolynomial(HCMK.support_coefficients(polysys6)[1])
    )
    PG6::ProjectiveGrid{Float64, UInt(3)} = ProjectiveGrid(
        gridPolySys6,
        UInt(3)
    )

    @test size(PG6.charts)[1] == 6
    @test isnothing(PG6.est_condition)

    projectiveGridHan!(PG6,UInt(1);maxDepth=UInt(10))

    # Test min depth increase
    oldMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    oldMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    projectiveIncreaseMinDepth!(PG6, UInt(5);scale=2.0)
    newMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    newMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    @test newMinDepthPG6==max(5,oldMinDepthPG6)
    @test newMaxDepthPG6==max(5,oldMaxDepthPG6)

    oldMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    oldMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    projectiveIncreaseDepth!(PG6, UInt(2);scale=2.0)
    newMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    newMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    @test newMinDepthPG6==oldMinDepthPG6+2
    @test newMaxDepthPG6==oldMaxDepthPG6+2

    oldMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    oldMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    projectiveIncreaseMinDepth!(PG6, UInt(8);scale=2.0,nodeFilter=(node)->depth(node)≤10)
    newMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    newMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    @test newMinDepthPG6==max(oldMinDepthPG6,8)
    @test newMaxDepthPG6≤10

    oldMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    oldMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    projectiveIncreaseDepth!(PG6, UInt(1);scale=2.0,nodeFilter=(node)->depth(node)≤10)
    newMinDepthPG6 = minimum(
        G -> minimum(depth, gridnodes(G)),
        PG6.charts
    )
    newMaxDepthPG6 = maximum(
        G -> maximum(depth, gridnodes(G)),
        PG6.charts
    )
    @test newMinDepthPG6==max(oldMinDepthPG6,9)
    @test newMaxDepthPG6≤10
end
