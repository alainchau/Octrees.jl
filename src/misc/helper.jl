using LinearAlgebra

function nearest_neighbors(X, p, r)
    neighbors = Int[]
    for i in 1:size(X, 2)
        if norm(X[:,i] - p) < r
            push!(neighbors, i)
        end
    end
    return neighbors
end
