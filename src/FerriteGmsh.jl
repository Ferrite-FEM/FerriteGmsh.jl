module FerriteGmsh

using Ferrite:
    Ferrite, Grid, Node, Vec
using Gmsh: Gmsh, gmsh

# Compat for Ferrite before v1.0
const FacetIndex = isdefined(Ferrite, :FacetIndex) ? Ferrite.FacetIndex : Ferrite.FaceIndex
const facets     = isdefined(Ferrite, :facets)     ? Ferrite.facets     : Ferrite.faces

if !isdefined(Ferrite, :SerendipityQuadraticHexahedron)
    const SerendipityQuadraticHexahedron = Ferrite.Cell{3,20,6}
    const QuadraticHexahedron = Ferrite.Cell{3,27,6}
else
    const SerendipityQuadraticHexahedron = Ferrite.SerendipityQuadraticHexahedron
    const QuadraticHexahedron = Ferrite.QuadraticHexahedron
end

"""
    gmshtoferritecell::Dict{String,DataType}

Map the gmsh element name -- as reported by `gmsh.model.mesh.getElementProperties` -- to
the `Ferrite` cell type describing the same geometry. Together with
[`gmshtoferriteperm`](@ref) this covers every cell type that `Ferrite` defines.
"""
const gmshtoferritecell = Dict{String,DataType}(
    "Line 2" => Ferrite.Line,
    "Line 3" => Ferrite.QuadraticLine,
    "Triangle 3" => Ferrite.Triangle,
    "Triangle 6" => Ferrite.QuadraticTriangle,
    "Quadrilateral 4" => Ferrite.Quadrilateral,
    "Quadrilateral 9" => Ferrite.QuadraticQuadrilateral,
    "Tetrahedron 4" => Ferrite.Tetrahedron,
    "Tetrahedron 10" => Ferrite.QuadraticTetrahedron,
    "Hexahedron 8" => Ferrite.Hexahedron,
    "Hexahedron 20" => SerendipityQuadraticHexahedron,
    "Hexahedron 27" => QuadraticHexahedron,
)

# Cells that only exist in newer Ferrite versions: `Wedge` was added in Ferrite 0.3.14,
# `Pyramid` and `SerendipityQuadraticQuadrilateral` in Ferrite 1.0, and `QuadraticWedge`
# after Ferrite 1.x. The node ordering of "Prism 18" matches `Lagrange{RefPrism, 2}`, so it
# needs no entry in `gmshtoferriteperm`.
for (gmshname, ferritename) in (("Quadrilateral 8", :SerendipityQuadraticQuadrilateral),
                                ("Prism 6", :Wedge),
                                ("Prism 18", :QuadraticWedge),
                                ("Pyramid 5", :Pyramid))
    isdefined(Ferrite, ferritename) && (gmshtoferritecell[gmshname] = getfield(Ferrite, ferritename))
end

"""
    gmshtoferriteperm::Dict{String,Tuple}

Local node permutation for the gmsh element types whose node ordering differs from the
corresponding `Ferrite` cell, such that

    ferrite_cell_nodes[i] == gmsh_element_nodes[gmshtoferriteperm[name][i]]

Element types absent from this `Dict` use the identity permutation.

Each entry is obtained by matching the gmsh local node coordinates -- the fifth return
value of `gmsh.model.mesh.getElementProperties` -- against
`Ferrite.reference_coordinates(Ferrite.geometric_interpolation(cell))`, up to the affine
map between the two reference domains. The vertices always come first in both codes, the
permutations only reorder within the vertex/edge/face/interior blocks. See
`test/test_cell_types.jl`, which pins every entry down by checking that a single element
of known shape integrates to its exact volume with a positive Jacobian.
"""
const gmshtoferriteperm = Dict{String,Tuple}(
    #  y
    #  ^
    #  |
    #  +--- > x
    #   \\
    #    z
    #
    # Gmsh numbers the two mid-edge nodes on the edges touching the apex the other way
    # round than Ferrite does.
    "Tetrahedron 10" => (1, 2, 3, 4, 5, 6, 7, 8, 10, 9),

    #     GMSH                  Ferrite
    # 4----14----3          4----11----3
    # |\         |\         |\         |\
    # |16        | 15       |20        | 19
    # 10 \       12 \       12 \       10 \
    # |   8----20+---7      |   8----15+---7
    # |   |      |   |      |   |      |   |
    # 1---+-9----2   |      1---+-9----2   |
    #  \ 18       \  19      \ 16       \  14
    #  11 |        13|       17 |        18|
    #    \|         \|         \|         \|
    #     5----17----6          5----13----6
    "Hexahedron 20" => (1, 2, 3, 4, 5, 6, 7, 8,
                        9, 12, 14, 10, 17, 19, 20, 18, 11, 13, 15, 16),

    # Same edge block as "Hexahedron 20"; the six face nodes 21-26 follow the face
    # ordering of the respective code (gmsh: -z, -y, -x, +x, +y, +z;
    # Ferrite: -z, -y, +x, +y, -x, +z), node 27 is the cell interior in both.
    "Hexahedron 27" => (1, 2, 3, 4, 5, 6, 7, 8,
                        9, 12, 14, 10, 17, 19, 20, 18, 11, 13, 15, 16,
                        21, 22, 24, 25, 23, 26, 27),

    # Gmsh walks the quadrilateral base counter-clockwise (1-2-3-4) while Ferrite orders
    # it lexicographically, (0,0)-(1,0)-(0,1)-(1,1), i.e. the cycle 1-2-4-3. Node 5 is
    # the apex in both.
    "Pyramid 5" => (1, 2, 4, 3, 5),
)

function _tocell(elementname::String)
    cell = get(gmshtoferritecell, elementname, nothing)
    cell === nothing && error("""
        unsupported gmsh element type "$elementname": Ferrite has no cell type with a \
        matching geometry. Supported gmsh element types are \
        $(join(sort!(collect(keys(gmshtoferritecell))), ", ")). Note that Ferrite only \
        provides first and second order cells, so higher order gmsh elements cannot be \
        converted -- re-mesh with `gmsh.model.mesh.setOrder(1)` or `setOrder(2)`.""")
    return cell
end

# Function barrier: `perm` comes out of an abstractly typed `Dict`, specializing here on
# its concrete `NTuple{N,Int}` type keeps the inner loop type stable.
function _tocells(::Type{CellType}, nodetags::Vector{Int64}, perm::NTuple{N,Int}) where {CellType,N}
    return [CellType(ntuple(j -> nodetags[i + perm[j] - 1], Val(N))) for i in 1:N:length(nodetags)]
end

"""
    cellorientation::Dict{DataType,Tuple{Tuple,Tuple}}

For each supported cell type, the pair `(simplex, flip)` used by [`reorient!`](@ref) to
detect and repair cells that Gmsh emitted with the opposite orientation.

`simplex` are `refdim + 1` local node indices spanning a *positively* oriented simplex of
the reference element. Taking the same determinant over the physical nodes therefore
yields the sign of `det(J)` at that corner of the cell.

`flip` is the relabelling induced by an orientation reversing symmetry `R` of the
reference element -- swapping the first two reference coordinates, or negating the
coordinate of a line. Relabelling a cell this way leaves it on exactly the same physical
region and only flips the sign of `det(J)`: for a nodal basis the geometric mapping turns
from `phi(xi)` into `phi(R(xi))`. Being induced by a reflection, `flip` is an involution.

`test/test_orientation.jl` rederives both entries of every cell type from
`Ferrite.reference_coordinates`.
"""
const cellorientation = Dict{DataType,Tuple{Tuple,Tuple}}()

for (gmshname, simplex, flip) in (
        ("Line 2", (1, 2), (2, 1)),
        ("Line 3", (1, 2), (2, 1, 3)),
        ("Triangle 3", (1, 2, 3), (2, 1, 3)),
        ("Triangle 6", (1, 2, 3), (2, 1, 3, 4, 6, 5)),
        ("Quadrilateral 4", (1, 2, 3), (1, 4, 3, 2)),
        ("Quadrilateral 8", (1, 2, 3), (1, 4, 3, 2, 8, 7, 6, 5)),
        ("Quadrilateral 9", (1, 2, 3), (1, 4, 3, 2, 8, 7, 6, 5, 9)),
        ("Tetrahedron 4", (1, 2, 3, 4), (1, 3, 2, 4)),
        ("Tetrahedron 10", (1, 2, 3, 4), (1, 3, 2, 4, 7, 6, 5, 8, 10, 9)),
        ("Hexahedron 8", (1, 2, 4, 5), (1, 4, 3, 2, 5, 8, 7, 6)),
        ("Hexahedron 20", (1, 2, 4, 5),
            (1, 4, 3, 2, 5, 8, 7, 6, 12, 11, 10, 9, 16, 15, 14, 13, 17, 20, 19, 18)),
        ("Hexahedron 27", (1, 2, 4, 5),
            (1, 4, 3, 2, 5, 8, 7, 6, 12, 11, 10, 9, 16, 15, 14, 13, 17, 20, 19, 18,
             21, 25, 24, 23, 22, 26, 27)),
        ("Prism 6", (1, 2, 3, 4), (1, 3, 2, 4, 6, 5)),
        ("Prism 18", (1, 2, 3, 4),
            (1, 3, 2, 4, 6, 5, 8, 7, 9, 10, 12, 11, 14, 13, 15, 17, 16, 18)),
        ("Pyramid 5", (1, 2, 3, 5), (1, 3, 2, 4, 5)),
    )
    haskey(gmshtoferritecell, gmshname) &&
        (cellorientation[gmshtoferritecell[gmshname]] = (simplex, flip))
end

_det(e::Tuple{Vec{1,T}}) where {T} = e[1][1]
_det(e::Tuple{Vec{2,T},Vec{2,T}}) where {T} = e[1][1] * e[2][2] - e[1][2] * e[2][1]
function _det(e::Tuple{Vec{3,T},Vec{3,T},Vec{3,T}}) where {T}
    a, b, c = e
    return a[1] * (b[2] * c[3] - b[3] * c[2]) -
           a[2] * (b[1] * c[3] - b[3] * c[1]) +
           a[3] * (b[1] * c[2] - b[2] * c[1])
end

# Function barriers: `simplex` and `flip` come out of an abstractly typed `Dict`.
function _simplexdet(cellnodes::NTuple, nodes, simplex::NTuple{M,Int}) where {M}
    x0 = nodes[cellnodes[simplex[1]]].x
    return _det(ntuple(k -> nodes[cellnodes[simplex[k + 1]]].x - x0, Val(M - 1)))
end
_flipnodes(cellnodes::NTuple{N,Int}, flip::NTuple{N,Int}) where {N} =
    ntuple(i -> cellnodes[flip[i]], Val(N))

"""
    reorient!(elements, nodes)

Relabel in place every cell of `elements` whose Jacobian determinant is negative, and
return `elements`.

Gmsh numbers the nodes of an element following the orientation of the entity the element
belongs to, so a surface whose normal points along `-z` yields clockwise elements while
Ferrite requires a positive Jacobian determinant. The relabelling keeps a cell on exactly
the same physical region, it only reverses its orientation, so the resulting grid
describes the same mesh. Cells that are already oriented correctly are left untouched.

This is applied automatically by [`togrid`](@ref); call it explicitly when assembling a
`Ferrite.Grid` from [`tonodes`](@ref) and [`toelements`](@ref) by hand, before computing
facet sets, since flipping a cell renumbers its local facets.
"""
function reorient!(elements::AbstractVector{<:Ferrite.AbstractCell}, nodes::AbstractVector{<:Node})
    nflipped = 0
    for (i, element) in pairs(elements)
        simplex, flip = get(cellorientation, typeof(element)) do
            error("cannot determine the orientation of a $(typeof(element)): unknown cell type")
        end
        if _simplexdet(element.nodes, nodes, simplex) < 0
            elements[i] = typeof(element)(_flipnodes(element.nodes, flip))
            nflipped += 1
        end
    end
    @debug "reoriented $nflipped of $(length(elements)) cells"
    return elements
end

function tonodes()
    nodeid, nodes = gmsh.model.mesh.getNodes()
    dim = Int64(gmsh.model.getDimension()) # Int64 otherwise julia crashes
    return [Node(Vec{dim}(nodes[i:i + (dim - 1)])) for i in 1:3:length(nodes)]
end

function toelements(dim::Int)
    elementtypes, elementtags, nodetags = gmsh.model.mesh.getElements(dim, -1)
    if isempty(elementtypes)
        error("could not find any elements with dimension $dim")
    end
    nodetags_all = convert(Vector{Vector{Int64}}, nodetags)
    if length(elementtypes) == 1
        elementname, _, _, _, _, _ = gmsh.model.mesh.getElementProperties(elementtypes[1])
        elements = _tocell(elementname)[]
    else
        elements = Ferrite.AbstractCell[]
    end

    for (eletypeidx,eletype) in enumerate(elementtypes)
        nodetags = nodetags_all[eletypeidx]
        elementname, dim, order, numnodes, localnodecoord, numprimarynodes = gmsh.model.mesh.getElementProperties(eletype)
        ferritecell = _tocell(elementname)
        perm = get(() -> ntuple(identity, numnodes), gmshtoferriteperm, elementname)
        elements_batch = _tocells(ferritecell, nodetags, perm)
        append!(elements,elements_batch)
    end

    return elements, reduce(vcat,convert(Vector{Vector{Int64}}, elementtags))
end

function toboundary(dim::Int)
    boundarydict = Dict{String,Vector}()
    boundaries = gmsh.model.getPhysicalGroups(dim) 
    for boundary in boundaries
        physicaltag = boundary[2]
        name = gmsh.model.getPhysicalName(dim, physicaltag)
        boundaryentities = gmsh.model.getEntitiesForPhysicalGroup(dim, physicaltag)
        boundaryconnectivity = Tuple[]
        for entity in boundaryentities
            boundarytypes, boundarytags, boundarynodetags = gmsh.model.mesh.getElements(dim, entity)
            boundarynodetags_all = convert(Vector{Vector{Int64}}, boundarynodetags)
            # A single entity can hold more than one element type, e.g. the surface of a
            # hybrid mesh carries both triangles and quadrilaterals, so collect them all.
            for (typeidx, boundarytype) in enumerate(boundarytypes)
                _, _, _, numnodes, _, _ = gmsh.model.mesh.getElementProperties(boundarytype)
                tags = boundarynodetags_all[typeidx]
                append!(boundaryconnectivity, [Tuple(tags[i:i + (numnodes - 1)]) for i in 1:numnodes:length(tags)])
            end
        end
        boundarydict[name] = boundaryconnectivity
    end 
    return boundarydict
end

function _add_to_facetsettuple!(facetsettuple::Set{FacetIndex}, boundaryfacet::Tuple, element_facets)
    for (eleidx, elefacets) in enumerate(element_facets)
        if any(issubset.(elefacets, (boundaryfacet,)))
            localfacet = findfirst(x -> issubset(x,boundaryfacet), elefacets) 
            push!(facetsettuple, FacetIndex(eleidx, localfacet))
        end
    end
    return facetsettuple
end

function tofacetsets(boundarydict::Dict{String,Vector}, elements::Vector{<:Ferrite.AbstractCell})
    element_facets = facets.(elements)
    facetsets = Dict{String,Set{FacetIndex}}()
    for (boundaryname, boundaryfacets) in boundarydict
        facetsettuple = Set{FacetIndex}()
        for boundaryfacet in boundaryfacets
            _add_to_facetsettuple!(facetsettuple, boundaryfacet, element_facets)
        end
        facetsets[boundaryname] = facetsettuple
    end
    return facetsets
end

function tocellsets(dim::Int, global_elementtags::Vector{Int})
    cellsets = Dict{String,Set{Int}}()
    element_to_cell_mapping = Dict(zip(global_elementtags, eachindex(global_elementtags)))
    physicalgroups = gmsh.model.getPhysicalGroups(dim)
    for (_, physicaltag) in physicalgroups 
        gmshname = gmsh.model.getPhysicalName(dim, physicaltag)
        isempty(gmshname) ? (name = "$physicaltag") : (name = gmshname)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim,physicaltag)
        cellsetelements = Set{Int}()
        for entity in entities
            _, elementtags, _= gmsh.model.mesh.getElements(dim, entity)
            elementtags = reduce(vcat,elementtags) |> x-> convert(Vector{Int},x)
            for ele in elementtags
                push!(cellsetelements, element_to_cell_mapping[ele])
            end
            cellsets[name] = cellsetelements
        end
    end
    return cellsets
end

"""
    togrid(filename::String; domain="")

Open the Gmsh file `filename` (ie a `.geo` or `.msh` file) and return the corresponding
`Ferrite.Grid`.
"""
function togrid(filename::String; domain="")
    # Check that file exists since Gmsh assumes we want to start a new model
    # if passing a non-existing path. In this function we need the model to exist.
    if !isfile(filename)
        # This is the error that is thrown by open("non-existing"),
        # error code 2 is "no such file or directory".
        throw(SystemError("opening file $(repr(filename))", 2))
    end
    should_finalize = Gmsh.initialize()
    gmsh.open(filename)
    fileextension = filename[findlast(isequal('.'), filename):end]
    dim = Int64(gmsh.model.getDimension()) # dont ask..

    if fileextension != ".msh"
        gmsh.model.mesh.generate(dim)
    end
    grid = togrid(; domain=domain)
    should_finalize && Gmsh.finalize()
    return grid
end

@deprecate saved_file_to_grid togrid

"""
    togrid(; domain="")

Generate a `Ferrite.Grid` from the current active/open model in the Gmsh library.
"""
function togrid(; domain="")
    dim = Int64(gmsh.model.getDimension())
    facedim = dim - 1
    saveall_flag = Bool(gmsh.option.getNumber("Mesh.SaveAll"))
    # set the save_all flag to one. hotfix #TODO for future
    if !saveall_flag
        gmsh.option.setNumber("Mesh.SaveAll",1)
    end
    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()
    nodes = tonodes()
    elements, gmsh_elementidx = toelements(dim)
    # Must happen before `tofacetsets` below, since flipping a cell renumbers its facets.
    reorient!(elements, nodes)
    cellsets = tocellsets(dim, gmsh_elementidx)

    if !isempty(domain)
        domaincellset = cellsets[domain]
        elements = elements[collect(domaincellset)]
    end

    boundarydict = toboundary(facedim)
    facetsets = tofacetsets(boundarydict, elements)
    # reset the save_all flag to the default value
    if !saveall_flag
        gmsh.option.setNumber("Mesh.SaveAll",0)
    end
    @static if isdefined(Ferrite, :FacetIndex)
        return Grid(elements, nodes, facetsets=facetsets, cellsets=cellsets)
    else # Compat for Ferrite before v1.0
        return Grid(elements, nodes, facesets=facetsets, cellsets=cellsets)
    end
end

export gmsh
export tonodes, toelements, toboundary, tofacetsets, tocellsets, togrid, reorient!

@deprecate tofacesets tofacetsets

end
