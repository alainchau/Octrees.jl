struct Octree
    X
    root        # Parent cell
    maxdepth
end

struct LocationalCode
    # TODO Should these bit strings be stored as symbols?
    x::String
    y::String
    z::String
end
LocationalCode() = LocationalCode("", "", "")

struct LeafData
    ids
    code
    depth
end
LeafData(;ids, code, depth) = LeafData(ids, code, depth)
