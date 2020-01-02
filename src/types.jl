struct Octree{T<:AbstractFloat}
    X::Array{T}
    root::Cell        # Parent cell
    maxdepth::Int
end

struct LocationalCode
    x::String
    y::String
    z::String
end
LocationalCode() = LocationalCode("", "", "")

struct LeafData
    ids::Array{Int}
    code::LocationalCode
    depth::Int
end

LeafData(; ids, code, depth) = LeafData(ids, code, depth)
