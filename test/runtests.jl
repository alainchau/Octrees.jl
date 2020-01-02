using Octrees
using Test
using LinearAlgebra
using Random

include("../src/misc/helper.jl") #import nearest neighbors

@testset "Octrees.jl" begin
    # getindex
    xloc = randstring("01", 10)
    yloc = randstring("01", 10)
    zloc = randstring("01", 10)
    code = LocationalCode(xloc, yloc, zloc)
    @test code[1:5] == LocationalCode(xloc[1:5], yloc[1:5], zloc[1:5])
    @test code[1:end] == LocationalCode(xloc, yloc, zloc)
    @test code[3:end] == LocationalCode(xloc[3:end], yloc[3:end], zloc[3:end])

    @test code * [1,0,0] == LocationalCode(xloc * "1", yloc * "0", zloc * "0")
    @test code * "100" == LocationalCode(xloc * "1", yloc * "0", zloc * "0")

    # Test locational code functionality
    N = 100
    X = randn(3, N)
    octree = Octree(X)
    for i in rand(1:N, 10)
        @test i in code_to_cell(octree, LocationalCode(octree, X[:, i])).data.ids
    end

    # Check if brute force checking for nearest neighbors is equivalent to
    # finding nearest neighbors via octree.
    for iteration in 1:10
        n = 100
        X = randn(3, n)
        p = [0, 0, 0.]
        r = 0.5
        octree = Octree(X, 0.05, 1)
        octree_neighbors = knn(octree, p, r)
        actual_neighbors = nearest_neighbors(X, p, r) |> Set
        @test octree_neighbors == actual_neighbors
    end
end
