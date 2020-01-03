# Octrees

[![Build Status](https://travis-ci.com/alainchau/Octrees.jl.svg?branch=master)](https://travis-ci.com/alainchau/Octrees.jl)
[![Codecov](https://codecov.io/gh/alainchau/Octrees.jl/branch/master/graph/badge.svg?token=a5St9etPNS)](https://codecov.io/gh/alainchau/Octrees.jl)

Use locational codes for random access of cells of an octree. In particular, this technique is useful for quickly finding the nearest neighbors of an arbitrary point in the octree.



# Installation
```julia
]add https://github.com/alainchau/Octrees.jl
```

# Example
```julia
julia> Using Octrees

# Generate artifical data and project onto xy-plane
julia> X = randn(3,100); X[3,:] .= 0

# Create octree
julia> octree = Octree(X);
Creating Octree with minimum side length δ = 0.02745
   
julia> using Plots
julia> include("src/misc/plotstuff.jl")
plot3d! (generic function with 2 methods)
   
julia> plot!(octree)
  
# Find nearest neighbors
julia> nn = knn(octree, [0,0,0], 1) |> collect;
julia> scatter!(X[1,nn], X[2,nn], markersize=3, color=:red)

# Draw circle to verify nearest neighbors
julia> ts = range(0.,2π,length=100)
julia> xs, ys = cos.(ts), sin.(ts)
julia> plot!(xs, ys, marker=0, fillcolor=:red, fillalpha=0.5, seriestype=:shape)
```

![alt text](https://github.com/alainchau/Octrees.jl/blob/master/src/misc/octree_example.png "Logo Title Text 1")


Relevant literature: http://ronaldperry.org/treeTraversalJGTWithCode.pdf
