# Plot stuff
import Plots: plot!, plot3d!

function plot!(octree::Octree)
    plt = plot(legend=nothing)
    for leaf in allleaves(octree.root)
       v = hcat(collect(vertices(leaf.boundary))...)
       plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]])
   end
   plt
end

function plot3d!(octree::Octree)
    plt = plot(legend=nothing)
    for leaf in allleaves(octree.root)
       v = hcat(collect(vertices(leaf.boundary))...)
       facets = [[1,2,4,3,1], [1,5,7,3], [2,6,8,4], [7,8], [5,6]]
       for facet in facets
           plot!(plt, v[1, facet], v[2, facet], v[3, facet])
       end
   end
   plt
end
