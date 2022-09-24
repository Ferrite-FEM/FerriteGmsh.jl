module FerriteGmsh
using Ferrite

import gmsh_jll
include(gmsh_jll.gmsh_api)
import .gmsh

"""
    FerriteGmsh.initialize(; finalize_atexit=true)

Wrapper around `gmsh.initialize` which make sure to only call it if `gmsh` is not already
initialized. Return `true` if `gmsh.initialize` was called, and `false` if `gmsh` was
already initialized.

If `finalize_atexit` is `true` a Julia exit hook is added, which calls
`gmsh.finalize()`.
"""
function initialize(; finalize_atexit=true)
    if Bool(gmsh.isInitialized())
        return false
    end
    gmsh.initialize()
    if finalize_atexit
        atexit(finalize)
    end
    return true
end

"""
    FerriteGmsh.finalize()

Wrapper around `gmsh.finalize` which make sure to only call it if `gmsh` is initialized.
"""
function finalize()
    if Bool(gmsh.isInitialized())
        gmsh.finalize()
    end
end

const gmshtoferritecell = Dict("Line 2" => Line,
                              "Line 3" => QuadraticLine,
                              "Triangle 3" => Triangle,
                              "Triangle 6" => QuadraticTriangle,
                              "Quadrilateral 4" => Quadrilateral,
                              "Quadrilateral 9" => QuadraticQuadrilateral,
                              "Tetrahedron 4" => Tetrahedron,
                              "Tetrahedron 10" => QuadraticTetrahedron,
                              "Hexahedron 8" => Hexahedron,
                              "Hexahedron 20" => Cell{3,20,6})

function translate_elements(original_elements)
    return original_elements
end

function translate_elements(original_elements::Vector{QuadraticTetrahedron})
    ferrite_elements = QuadraticTetrahedron[]
    for original_ele in original_elements
        push!(ferrite_elements,QuadraticTetrahedron((original_ele.nodes[1], 
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

function translate_elements(original_elements::Vector{Cell{3,20,6}})
    ferrite_elements = Cell{3,20,6}[]
    for original_ele in original_elements
        push!(ferrite_elements,Cell{3,20,6}((original_ele.nodes[1], 
                                             original_ele.nodes[2], 
                                             original_ele.nodes[3], 
                                             original_ele.nodes[4],
                                             original_ele.nodes[5],  
                                             original_ele.nodes[6], 
                                             original_ele.nodes[7], 
                                             original_ele.nodes[8],
                                             original_ele.nodes[9],
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

function tofacesets(boundarydict::Dict{String,Vector}, elements::Vector{<:Ferrite.AbstractCell})
    faces = Ferrite.faces.(elements)
    facesets = Dict{String,Set{FaceIndex}}()
    for (boundaryname, boundaryfaces) in boundarydict
        facesettuple = Set{FaceIndex}()
        for boundaryface in boundaryfaces
            for (eleidx, elefaces) in enumerate(faces)
                if any(issubset.(elefaces, (boundaryface,)))
                    localface = findfirst(x -> issubset(x,boundaryface), elefaces) 
                    push!(facesettuple, FaceIndex(eleidx, localface))
                end
            end
        end
        facesets[boundaryname] = facesettuple
    end
    return facesets
end

function tocellsets(dim::Int, global_elementtags::Vector{Int})
    cellsets = Dict{String,Set{Int}}()
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
                push!(cellsetelements, findfirst(isequal(ele), global_elementtags))
            end
            cellsets[name] = cellsetelements
        end
    end
    return cellsets
end

function saved_file_to_grid(filename::String; domain="")
    # Check that file exists since Gmsh assumes we want to start a new model
    # if passing a non-existing path. In this function we need the model to exist.
    if !isfile(filename)
        # This is the error that is thrown by open("non-existing"),
        # error code 2 is "no such file or directory".
        throw(SystemError("opening file $(repr(filename))", 2))
    end
    should_finalize = initialize()
    gmsh.open(filename)
    fileextension = filename[findlast(isequal('.'), filename):end]
    dim = Int64(gmsh.model.getDimension()) # dont ask..
    facedim = dim - 1

    if fileextension != ".msh"
        gmsh.model.mesh.generate(dim)
    end
 
    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()
    nodes = tonodes()
    elements, gmsh_elementidx = toelements(dim) 
    cellsets = tocellsets(dim, gmsh_elementidx)

    if !isempty(domain)
        domaincellset = cellsets[domain]
        elements = elements[collect(domaincellset)]
    end

    boundarydict = toboundary(facedim)
    facesets = tofacesets(boundarydict, elements)

    should_finalize && finalize()

    return Grid(elements, nodes, facesets=facesets, cellsets=cellsets)
end

export gmsh
export tonodes, toelements, toboundary, tofacesets, tocellsets, saved_file_to_grid

end
