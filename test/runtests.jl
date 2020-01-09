using Octrees
using Test
using LinearAlgebra
using Random
using Plots
include("../src/misc/plotstuff.jl")
include("../src/misc/helper.jl")    # import nearest neighbors

@testset "Octrees.jl" begin
    # getindex
    xloc = randstring("01", 10)
    yloc = randstring("01", 10)
    zloc = randstring("01", 10)

    # constructor
    code = LocationalCode(xloc, yloc, zloc)
    @test code[1:5] == LocationalCode(xloc[1:5], yloc[1:5], zloc[1:5])
    @test code[1:end] == LocationalCode(xloc, yloc, zloc)
    @test code[3:end] == LocationalCode(xloc[3:end], yloc[3:end], zloc[3:end])

    # ops
    @test code * [1,0,0] == LocationalCode(xloc * "1", yloc * "0", zloc * "0")
    @test code * "100" == LocationalCode(xloc * "1", yloc * "0", zloc * "0")

    # Test locational code functionality
    for N in [100, 150, 200]
        X = randn(3, N)
        octree = Octree(X)
        for i in rand(1:N, 10)
            @test i in code_to_cell(octree, LocationalCode(octree, X[:, i])).data.ids
        end
    end
    
    # Check if brute force checking for nearest neighbors is equivalent to
    # finding nearest neighbors via octree.
    for iteration in 1:10
        n = 100
        X = randn(3, n)
        p = X[:, 1]
        r = 0.5
        octree = Octree(X)
        octree_neighbors = knn(octree, p, r)
        actual_neighbors = nearest_neighbors(X, p, r) |> Set
        @test octree_neighbors == actual_neighbors
    end

    # Sample from sphere
    for iteration in 1:10
        n = 300
        X = randn(3, n)
        X ./= mapslices(norm, X, dims=1)
        p = X[:, 1]
        r = 0.2
        octree = Octree(X)
        octree_neighbors = knn(octree, p, r)
        actual_neighbors = nearest_neighbors(X, p, r) |> Set
        @test octree_neighbors == actual_neighbors
    end

    # Plot octree of data projected onto xy-plane.
    X = randn(3,100)
    X[3,:] .= 0
    plt = plot!(Octree(X))
    savefig(plt, "testoutput.png")
end
