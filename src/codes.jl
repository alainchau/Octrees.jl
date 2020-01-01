
function *(loc::LocationalCode, bits::String)
    xloc = loc.x * bits[1]
    yloc = loc.y * bits[2]
    zloc = loc.z * bits[3]
    return LocationalCode(xloc, yloc, zloc)
end

*(loc::LocationalCode, bits::Vector{Int}) = *(loc, join(string.(bits)))
*(loc::LocationalCode, bits::Tuple{Int, Int, Int}) = *(loc, join(string.(bits)))

function LocationalCode(octree::Octree, p)
    binsize = 2^octree.maxdepth
    octree_size = octree.root.boundary.widths[1]

    # Convert to code
    xloc = floor(Int, (p[1] - octree.root.boundary.origin[1]) / octree_size * binsize)
    yloc = floor(Int, (p[2] - octree.root.boundary.origin[2]) / octree_size * binsize)
    zloc = floor(Int, (p[3] - octree.root.boundary.origin[3]) / octree_size * binsize)

    # Convert to strings
    xloc = bitstring(xloc)[end-octree.maxdepth+1:end]
    yloc = bitstring(yloc)[end-octree.maxdepth+1:end]
    zloc = bitstring(zloc)[end-octree.maxdepth+1:end]
    LocationalCode(xloc, yloc, zloc)
end

function Base.getindex(code::LocationalCode, inds...)
    return LocationalCode(code.x[inds...], code.y[inds...], code.z[inds...])
end

Base.lastindex(code::LocationalCode) = Base.lastindex(code.x)

Base.iterate(code::LocationalCode) = ((parse(Int, code.x[1]),
                                       parse(Int, code.y[1]),
                                       parse(Int, code.z[1])),
                                       2)
function Base.iterate(code::LocationalCode, state)
    if state <= length(code.x)
        return ((parse(Int, code.x[state]),
                 parse(Int, code.y[state]),
                 parse(Int, code.z[state])),
                 state+1)
    else
        return nothing
    end
end

"""Check if code is blank (corresponding to root cell)."""
isroot(code::LocationalCode) = (code == LocationalCode("", "", ""))

function code_to_cell(octree::Octree, code::LocationalCode)
    cell = octree.root

    # If cell is root, then return
    isroot(code) && return cell

    for (i, (xloc, yloc, zloc)) in enumerate(code)
        i > octree.maxdepth && break
        cell = cell[xloc+1, yloc+1, zloc+1]
        isnothing(cell.children) && break
    end
    return cell
end

toint(code::LocationalCode) = parse(Int, code.x, base=2)
