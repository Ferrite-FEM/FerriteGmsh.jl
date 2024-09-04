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

const gmshtoferritecell = Dict("Line 2" => Ferrite.Line,
                              "Line 3" => Ferrite.QuadraticLine,
                              "Triangle 3" => Ferrite.Triangle,
                              "Triangle 6" => Ferrite.QuadraticTriangle,
                              "Quadrilateral 4" => Ferrite.Quadrilateral,
                              "Quadrilateral 9" => Ferrite.QuadraticQuadrilateral,
                              "Tetrahedron 4" => Ferrite.Tetrahedron,
                              "Tetrahedron 10" => Ferrite.QuadraticTetrahedron,
                              "Hexahedron 8" => Ferrite.Hexahedron,
                              "Hexahedron 20" => SerendipityQuadraticHexahedron,
                              "Hexahedron 27"=> QuadraticHexahedron)

function translate_elements(original_elements)
    return original_elements
end

function translate_elements(original_elements::Vector{Ferrite.QuadraticTetrahedron})
    ferrite_elements = Ferrite.QuadraticTetrahedron[]
    for original_ele in original_elements
        push!(ferrite_elements,Ferrite.QuadraticTetrahedron((original_ele.nodes[1], 
                                                     original_ele.nodes[2], 
                                                     original_ele.nodes[3], 
                                                     original_ele.nodes[4],
                                                     original_ele.nodes[5],  
                                                     original_ele.nodes[6], 
                                                     original_ele.nodes[7], 
                                                     original_ele.nodes[8],
                                                     original_ele.nodes[10],
                                                     original_ele.nodes[9])),)
    end
    return ferrite_elements
end

#
#  y
#
#  ^
#  |
#  |
#  +--- > x
#  \\
#   \\
#    \\
#     z
#
#     GMSH
# 4----14----3
# |\         |\
# |16        | 15
# 10 \       12 \
# |   8----20+---7
# |   |      |   |
# 1---+-9----2   |
#  \ 18       \  19
#  11 |        13|
#    \|         \|
#     5----17----6
#
#    Ferrite
# 4----11----3
# |\         |\
# |20        | 19
# 12 \       10 \
# |   8----15+---7
# |   |      |   |
# 1---+-9----2   |
#  \ 16       \  14
#  17 |        18|
#    \|         \|
#     5----13----6
#
function translate_elements(original_elements::Vector{SerendipityQuadraticHexahedron})
    ferrite_elements = SerendipityQuadraticHexahedron[]
    for original_ele in original_elements
        push!(ferrite_elements, SerendipityQuadraticHexahedron((
                                             original_ele.nodes[1], 
                                             original_ele.nodes[2], 
                                             original_ele.nodes[3], 
                                             original_ele.nodes[4],
                                             original_ele.nodes[5],  
                                             original_ele.nodes[6], 
                                             original_ele.nodes[7], 
                                             original_ele.nodes[8],
                                             original_ele.nodes[9],  # edges
                                             original_ele.nodes[12],
                                             original_ele.nodes[14],
                                             original_ele.nodes[10],
                                             original_ele.nodes[17],
                                             original_ele.nodes[19],
                                             original_ele.nodes[20],
                                             original_ele.nodes[18],
                                             original_ele.nodes[11],
                                             original_ele.nodes[13],
                                             original_ele.nodes[15],
                                             original_ele.nodes[16])),)
    end
    return ferrite_elements
end


#  y
#
#  ^
#  |
#  |
#  +--- > x
#  \\
#   \\
#    \\
#     z
#
#     GMSH
# 4----14----3
# |\         |\
# |16    25  | 15
# 10 \ 21    12 \
# |   8----20+---7
# |23 |  27  | 24|
# 1---+-9----2   |
#  \ 18    26 \  19
#  11 |  22    13|
#    \|         \|
#     5----17----6
#
#    Ferrite
# 4----11----3
# |\         |\
# |20    24  | 19
# 12 \ 21    10 \
# |   8----15+---7
# |25 |  27  | 23|
# 1---+-9----2   |
#  \ 16     26\  14
#  17 |  22    18|
#    \|         \|
#     5----13----6
#
function translate_elements(original_elements::Vector{QuadraticHexahedron})
    ferrite_elements = QuadraticHexahedron[]
    for original_ele in original_elements
        push!(ferrite_elements,QuadraticHexahedron((
                                             original_ele.nodes[1], 
                                             original_ele.nodes[2], 
                                             original_ele.nodes[3], 
                                             original_ele.nodes[4],
                                             original_ele.nodes[5],  
                                             original_ele.nodes[6], 
                                             original_ele.nodes[7], 
                                             original_ele.nodes[8],
                                             original_ele.nodes[9], # edge
                                             original_ele.nodes[12],
                                             original_ele.nodes[14],
                                             original_ele.nodes[10],
                                             original_ele.nodes[17],
                                             original_ele.nodes[19],
                                             original_ele.nodes[20],
                                             original_ele.nodes[18],
                                             original_ele.nodes[11],
                                             original_ele.nodes[13],
                                             original_ele.nodes[15],
                                             original_ele.nodes[16],
                                             original_ele.nodes[21], # face
                                             original_ele.nodes[22],
                                             original_ele.nodes[24],
                                             original_ele.nodes[25],
                                             original_ele.nodes[23],
                                             original_ele.nodes[26],
                                             original_ele.nodes[27],
                                             )),)
    end
    return ferrite_elements
end 

function tonodes()
    nodeid, nodes = gmsh.model.mesh.getNodes()
    dim = Int64(gmsh.model.getDimension()) # Int64 otherwise julia crashes
    return [Node(Vec{dim}(nodes[i:i + (dim - 1)])) for i in 1:3:length(nodes)], Int64.(nodeid)
end

function toelements(dim::Int)
    elementtypes, elementtags, nodetags = gmsh.model.mesh.getElements(dim, -1)
    if isempty(elementtypes)
        error("could not find any elements with dimension $dim")
    end
    nodetags_all = convert(Vector{Vector{Int64}}, nodetags)
    if length(elementtypes) == 1
        elementname, _, _, _, _, _ = gmsh.model.mesh.getElementProperties(elementtypes[1])
        elements = gmshtoferritecell[elementname][]
    else
        elements = Ferrite.AbstractCell[]
    end

    for (eletypeidx,eletype) in enumerate(elementtypes)
        nodetags = nodetags_all[eletypeidx]
        elementname, dim, order, numnodes, localnodecoord, numprimarynodes = gmsh.model.mesh.getElementProperties(eletype) 
        ferritecell = gmshtoferritecell[elementname]
        elements_gmsh = [ferritecell(Tuple(nodetags[i:i + (numnodes - 1)])) for i in 1:numnodes:length(nodetags)]
        elements_batch = translate_elements(elements_gmsh)
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
            _, _, _, numnodes, _, _ = gmsh.model.mesh.getElementProperties(boundarytypes[1]) 
            boundarynodetags = convert(Vector{Vector{Int64}}, boundarynodetags)[1]
            append!(boundaryconnectivity, [Tuple(boundarynodetags[i:i + (numnodes - 1)]) for i in 1:numnodes:length(boundarynodetags)])
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

function todimEntitysets(dim::Int, global_entity_tags::Vector{Int})
    entityset = Dict{String,Set{Int}}()
    gmsh_to_ferrite_mapping = Dict(zip(global_entity_tags, eachindex(global_entity_tags)))
    physicalgroups = gmsh.model.getPhysicalGroups(dim)
    for (_, physicaltag) in physicalgroups 
        gmshname = gmsh.model.getPhysicalName(dim, physicaltag)
        isempty(gmshname) ? (name = "$physicaltag") : (name = gmshname)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim,physicaltag)
        ferrite_entities = Set{Int}()
        for entity in entities
            _, elementtags, _= gmsh.model.mesh.getElements(dim, entity)
            elementtags = reduce(vcat,elementtags) |> x-> convert(Vector{Int},x)
            for ele in elementtags
                push!(ferrite_entities, gmsh_to_ferrite_mapping[ele])
            end
            entityset[name] = ferrite_entities
        end
    end
    return entityset
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
    nodes, gmsh_nodeidx = tonodes()
    nodesets = todimEntitysets(0, gmsh_nodeidx)
    elements, gmsh_elementidx = toelements(dim) 
    cellsets = todimEntitysets(dim, gmsh_elementidx)

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
        return Grid(elements, nodes, facetsets=facetsets, cellsets=cellsets, nodesets=nodesets)
    else # Compat for Ferrite before v1.0
        return Grid(elements, nodes, facesets=facetsets, cellsets=cellsets, nodesets=nodesets)
    end
end

export gmsh
export tonodes, toelements, toboundary, tofacetsets, tocellsets, togrid

@deprecate tofacesets tofacetsets

end
