using Test
include("Grid.jl")


@testset "Grid" failfast=true begin
    grid = Grid(Float64, 2, Nothing)

    @test typeof(est_condition(grid)) == Float64
    @test dim(grid) == 2
end