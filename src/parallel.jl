# Read this and parallelize this stuff
# https://techytok.com/multiprocessing-in-julia-module/

function Base.display(cell::Cell)
    print("origin = $(round.(cell.boundary.origin, digits=3)).")
    print("\twidth =  $(round.(cell.boundary.widths[1], digits=3)).")
    print("\tisleaf = $(isnothing(cell.children))\n")
    println("   ids = $(cell.data.ids)")
end

function Base.display(octree::Octree)
    print("origin = $(round.(octree.root.boundary.origin, digits=3)).")
    print("\twidth =  $(round.(octree.root.boundary.widths[1], digits=3)).")
    println("\tisleaf = $(isnothing(octree.root.children)). depth=$(octree.maxdepth)")
end

abstract type MyRefinery <: AbstractRefinery end

struct SeqRefinery <: MyRefinery
    tolerance::Float64
    K::Int
    X::AbstractArray
end

struct ParRefinery <: MyRefinery
    tolerance::Float64
    K::Int
    X::AbstractArray
end

function in(x::Array, r::HyperRectangle{N,T}) where {N,T<:AbstractFloat}
    for i in eachindex(x)
        x[i] < r.origin[i] && return false
        x[i] > r.origin[i] + r.widths[i] && return false
    end
    return true
end

function needs_refinement(r::MyRefinery, cell)
    maximum(cell.boundary.widths) <= r.tolerance && return false
    length(cell.data.ids) > r.K
end

function refine_data(r::SeqRefinery, cell::Cell, indices)
    newpoints = Int32[]    
    for i in 1:length(cell.data.ids)
        if r.X[:, cell.data.ids[i]] in child_boundary(cell, indices)
            push!(newpoints, cell.data.ids[i])
        end
    end        
    return LeafData(ids = newpoints,
                    code = cell.data.code * (indices .- 1),
                    depth = cell.data.depth + 1)
end

function refine_data(r::ParRefinery, cell::Cell, indices)
    newpoints = Int32[]    
    @sync @distributed for i in 1:length(cell.data.ids)
        if r.X[:, cell.data.ids[i]] in child_boundary(cell, indices)
            push!(newpoints, cell.data.ids[i])
        end
    end        
    return LeafData(ids = newpoints,
                    code = cell.data.code * (indices .- 1),
                    depth = cell.data.depth + 1)
end


"""
    Octree(X::Array{S}, δ::T; K::Int=4) where {S<:M, T<:M} where M <: AbstractFloat

Create a quadtree, octree, or higher D-dimensional analog for the given dataset X with
dimensions D x N. The minimum number of points per cell is given by K and the minimum
side length of each cell is given by δ.
"""
function Octree(X::Array{S}, δ::T, K::Int, execution=:seq) where {S<:M,T<:M} where {M<:AbstractFloat}
    if execution == :seq
        r = SeqRefinery(δ, K, X)
    else
        r = ParRefinery(δ, K, X)
    end
    
    # Subtract small value to avoid corner cases
    origin = SVector(minimum(X, dims = 2) .- δ/2...)
    # Make a square bounding box
    widths = SVector(fill(
        maximum(maximum(X, dims = 2) - minimum(X, dims = 2)) .+ δ/2,
        size(X, 1),
    )...)
    root = Cell(origin, widths, LeafData(collect(1:size(X, 2)), LocationalCode(), Int16(0)))
    print("Making tree..")
    adaptivesampling!(root, r)
    println(" finished.")
    
    # Find max depth
    print("Finding max depth..")
    maxdepth = 0
    for leaf in allleaves(root)
        print(leaf.data.depth, " ")
        maxdepth = max(maxdepth, leaf.data.depth)
    end
    # maxdepth = maximum(map(x -> x.data.depth, allleaves(root)))
    println(" finished.")
    
    return Octree(X, root, maxdepth, execution)
end

function Octree(X::Array{S}) where S<:AbstractFloat
    """
    Use a heuristic to find appropriate values of δ and K.
    """
    δ = Inf
    for i in rand(1:size(X, 2), floor(Int, size(X, 2) / 2))
        for j in rand(1:size(X, 2), floor(Int, size(X, 2) / 2))
            i == j && continue
            δ = min(δ, norm(X[:, i] - X[:, j]))
        end
    end
    K = ceil(Int, log10(size(X, 2)))
    println("Creating Octree with minimum side length δ = $(round(δ, digits=5)) and K=$K.")
    return Octree(X, δ, K)
end

const NEIGHBOR_INDICES = (
    (-1, -1, -1),
    (-1, -1, 0),
    (-1, -1, 1),
    (-1, 0, -1),
    (-1, 0, 0),
    (-1, 0, 1),
    (-1, 1, -1),
    (-1, 1, 0),
    (-1, 1, 1),
    (0, -1, -1),
    (0, -1, 0),
    (0, -1, 1),
    (0, 0, -1),
    (0, 0, 0),
    (0, 0, 1),
    (0, 1, -1),
    (0, 1, 0),
    (0, 1, 1),
    (1, -1, -1),
    (1, -1, 0),
    (1, -1, 1),
    (1, 0, -1),
    (1, 0, 0),
    (1, 0, 1),
    (1, 1, -1),
    (1, 1, 0),
    (1, 1, 1),
)

"""
    knn(octree::Octree, p, r)

Find neighbors with fixed radius r of a given point p using an octree.
"""
function knn(octree::Octree, p::Vector, r::T) where T<:AbstractFloat
    ℓ = max(floor(Int, log(2, octree.root.boundary.widths[1] / r)), 1)
    code = LocationalCode(octree::Octree, p)
    neighbors = Set{Int}()

    # TODO Do i have to duplicate this? Maybe use macro?
    if octree.execution == :par
        @sync @distributed for (δ1, δ2, δ3) in NEIGHBOR_INDICES
            xloc = clamp(parse(Int, code.x[1:ℓ], base = 2) + δ1, 0, 2^(ℓ) - 1)
            yloc = clamp(parse(Int, code.y[1:ℓ], base = 2) + δ2, 0, 2^(ℓ) - 1)
            zloc = clamp(parse(Int, code.z[1:ℓ], base = 2) + δ3, 0, 2^(ℓ) - 1)
            xloc = bitstring(xloc)[end-ℓ+1:end]
            yloc = bitstring(yloc)[end-ℓ+1:end]
            zloc = bitstring(zloc)[end-ℓ+1:end]
            for i in code_to_cell(octree, LocationalCode(xloc, yloc, zloc)).data.ids
                if norm(octree.X[:, i] - p) < r
                    push!(neighbors, i)
                end
            end
        end
    else
        for (δ1, δ2, δ3) in NEIGHBOR_INDICES
            xloc = clamp(parse(Int, code.x[1:ℓ], base = 2) + δ1, 0, 2^(ℓ) - 1)
            yloc = clamp(parse(Int, code.y[1:ℓ], base = 2) + δ2, 0, 2^(ℓ) - 1)
            zloc = clamp(parse(Int, code.z[1:ℓ], base = 2) + δ3, 0, 2^(ℓ) - 1)
            xloc = bitstring(xloc)[end-ℓ+1:end]
            yloc = bitstring(yloc)[end-ℓ+1:end]
            zloc = bitstring(zloc)[end-ℓ+1:end]
            for i in code_to_cell(octree, LocationalCode(xloc, yloc, zloc)).data.ids
                if norm(octree.X[:, i] - p) < r
                    push!(neighbors, i)
                end
            end
        end
    end

    return neighbors
end
