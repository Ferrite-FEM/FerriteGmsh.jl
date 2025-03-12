function fix_node_order!(grid::Grid{2})
    cellids_from_type = Dict{Type, Vector{Int}}()
    for (i, cell) in enumerate(Ferrite.getcells(grid))
        ids = get!(cellids_from_type, typeof(cell), Int[])
        push!(ids, i)
    end

    for ids in values(cellids_from_type)
        _fix_node_order!(grid, ids)
    end
    return grid
end

function is_reversed(cell_coords::AbstractVector{<:Vec{2}})
    v1 = cell_coords[2] - cell_coords[1]
    v2 = cell_coords[3] - cell_coords[1]
    return last(Ferrite.:×(v1, v2)) < 0
end

function _fix_node_order!(grid::Grid{2}, ids::Vector{Int})
    x = Ferrite.getcoordinates(grid, first(ids))
    vertex_mapping = reverse_vertex_mapping(Ferrite.getcells(grid, first(ids)))
    facet_mapping = reverse_facet_mapping(Ferrite.getcells(grid, first(ids)))
    reversed_cells = Set{Int}()
    for i in ids
        Ferrite.getcoordinates!(x, grid, i)
        if is_reversed(x)
            grid.cells[i] = reverse_cell(grid.cells[i])
            push!(reversed_cells, i)
        end
    end
    if !isempty(reversed_cells)
        for set in values(grid.vertexsets)
            update_set!(set, reversed_cells, vertex_mapping)
        end
        for set in values(grid.facetsets)
            update_set!(set, reversed_cells, facet_mapping)
        end
    end
end

reverse_vertex_mapping(cell::Ferrite.AbstractCell) = reverse_vertex_mapping(Ferrite.getrefshape(cell))
reverse_facet_mapping(cell::Ferrite.AbstractCell) = reverse_facet_mapping(Ferrite.getrefshape(cell))

reverse_cell(cell::Ferrite.Triangle) = Ferrite.Triangle(reverse(cell.nodes))
function reverse_cell(cell::Ferrite.QuadraticTriangle)
    n = cell.nodes
    return Ferrite.QuadraticTriangle((n[3], n[2], n[1], n[5], n[4], n[6]))
end
reverse_vertex_mapping(::Type{Ferrite.RefTriangle}) = (3, 2, 1)
reverse_facet_mapping(::Type{Ferrite.RefTriangle}) = (2, 1, 3)

reverse_cell(cell::Ferrite.Quadrilateral) = Ferrite.Quadrilateral(reverse(cell.nodes))

function reverse_cell(cell::Ferrite.QuadraticQuadrilateral)
    n = cell.nodes
    return Ferrite.QuadraticQuadrilateral((n[4], n[3], n[2], n[1], n[7], n[6], n[5], n[8], n[9]))
end

function reverse_cell(cell::Ferrite.SerendipityQuadraticQuadrilateral)
    n = cell.nodes
    return Ferrite.SerendipityQuadraticQuadrilateral((n[4], n[3], n[2], n[1], n[7], n[6], n[5], n[8]))
end
reverse_vertex_mapping(::Type{Ferrite.RefQuadrilateral}) = (4, 3, 2, 1)
reverse_facet_mapping(::Type{Ferrite.RefQuadrilateral}) = (3, 2, 1, 4)

function update_set!(set::AbstractSet{BI}, reversed_cells, mapping) where {BI}
    # Note, this doesn't preserve ordering...
    new_items = Set{BI}()
    for item in set
        (cellnr, entitynr) = item
        if cellnr ∈ reversed_cells
            delete!(set, item)
            push!(new_items, BI(cellnr, mapping[entitynr]))
        end
    end
    union!(set, new_items)
    return set
end
