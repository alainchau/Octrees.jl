include("../src/misc/helper.jl")
using Octrees

X = randn(3, 5000)
octree = Octree(X)

nearest_neighbors(X, [0,0,0], 1.)
knn(octree, [0,0,0], 1.)

@time nearest_neighbors(X, [0,0,0], 1.)
@time knn(octree, [0,0,0], 1.)
