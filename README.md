# Octrees

[![Build Status](https://travis-ci.com/alainchau/Octrees.jl.svg?branch=master)](https://travis-ci.com/alainchau/Octrees.jl)
[![Codecov](https://codecov.io/gh/alainchau/Octrees.jl/branch/master/graph/badge.svg?token=a5St9etPNS)](https://codecov.io/gh/alainchau/Octrees.jl)

Use locational codes for random access of cells of an octree. In particular, this technique is useful for quickly finding the nearest neighbors of an arbitrary point in the octree.



# Installation
```julia
]add https://github.com/alainchau/Octrees.jl
```

# Example
## Plot
```julia
julia> Using Octrees

# Generate artifical data and project onto xy-plane
julia> X = randn(3,100); X[3,:] .= 0

# Create octree
julia> octree = Octree(X)
Creating Octree with minimum side length δ = 0.02745
origin = [-2.641, -2.196, -0.0].    width =  5.647. isleaf = false
   ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
   
julia> using Plots
julia> include("src/misc/plotstuff.jl")
plot3d! (generic function with 2 methods)
   
julia> plot!(octree)
  
# Find nearest neighbors
julia> nn = knn(octree, [0,0,0], 1) |> collect;
julia> scatter!(X[1,nn], X[2,nn], markersize=3, color=:red)
julia> scatter!(X[1,nn], X[2,nn], markersize=3, color=:red)

# Draw circle to verify nearest neighbors
julia> ts = range(0.,2π,length=100)
julia> xs, ys = cos.(ts), sin.(ts)
julia> plot!(xs, ys, marker=0, fillcolor=:red, fillalpha=0.5, seriestype=:shape)
```

![alt text](https://github.com/alainchau/Octrees.jl/blob/master/src/misc/octree_example.png "Logo Title Text 1")


Relevant literature: http://ronaldperry.org/treeTraversalJGTWithCode.pdf
