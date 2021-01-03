module JuAFEMGmsh

using gmsh
using JuAFEM

const gmshtojuafemcell = Dict("Line 2" => Line,
                              "Line 3" => QuadraticLine,
                              "Triangle 3" => Triangle,
                              "Triangle 6" => QuadraticTriangle,
                              "Quadrilateral 4" => Quadrilateral,
                              "Quadrilateral 9" => QuadraticQuadrilateral,
                              "Tetrahedron 4" => Tetrahedron,
                              "Tetrahedron 10" => QuadraticTetrahedron,
                              "Hexahedron 8" => Hexahedron,
                              "Hexahedron 20" => QuadraticHexahedron)

function tonodes()
    nodeid, nodes = gmsh.model.mesh.getNodes()
    dim = Int64(gmsh.model.getDimension()) # Int64 otherwise julia crashes
    return [Node(Vec{dim}(nodes[i:i + (dim - 1)])) for i in 1:3:length(nodes)]
end

function toelements(dim::Int)
    elementtypes, elementtags, nodetags = gmsh.model.mesh.getElements(dim, -1)
    @assert length(elementtypes) == 1 "only one element type per mesh is supported"
    elementname, dim, order, numnodes, localnodecoord, numprimarynodes = gmsh.model.mesh.getElementProperties(elementtypes[1]) 
    nodetags = convert(Array{Array{Int64,1},1}, nodetags)[1]
    juafemcell = gmshtojuafemcell[elementname]
    elements = [juafemcell(Tuple(nodetags[i:i + (numnodes - 1)])) for i in 1:numnodes:length(nodetags)]
    return elements, convert(Vector{Vector{Int64}}, elementtags)[1]
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

function tofacesets(boundarydict::Dict{String,Vector}, elements::Vector{<:JuAFEM.AbstractCell})
    faces = JuAFEM.faces.(elements)
    facesets = Dict{String,Set{Tuple{Int,Int}}}()
    for (boundaryname, boundaryfaces) in boundarydict
        facesettuple = Set{Tuple{Int,Int}}()
        for boundaryface in boundaryfaces
            for (eleidx, elefaces) in enumerate(faces)
                if boundaryface in elefaces
                    localface = findfirst(x -> x == boundaryface, elefaces) 
                    push!(facesettuple, (eleidx, localface))
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
        cellsetentities = gmsh.model.getEntitiesForPhysicalGroup(dim, physicaltag)
        cellsetelements = Set{Int}()
        for entity in cellsetentities
            _, elementtags, _ = gmsh.model.mesh.getElements(dim, entity)
            elementtags = convert(Vector{Vector{Int64}}, elementtags)[1]
            for elementtag in elementtags
                push!(cellsetelements, findfirst(isequal(elementtag), global_elementtags))
            end
        end
        cellsets[name] = cellsetelements
    end
    return cellsets
end

function saved_file_to_grid(filename::String; domain="")
    gmsh.initialize()
    gmsh.open(filename)
    fileextension = filename[findlast(isequal('.'), filename):end]
    dim = Int64(gmsh.model.getDimension()) # dont ask..
    facedim = dim - 1

    if fileextension != ".msh"
        gmsh.model.mesh.generate(dim)
    end    

    nodes = tonodes()
    elements, gmsh_elementidx = toelements(dim) 
    boundarydict = toboundary(facedim)
    facesets = tofacesets(boundarydict, elements)
    cellsets = tocellsets(dim, gmsh_elementidx)
    gmsh.finalize()

    if !isempty(domain)
        domaincellset = cellsets[domain]
        elements = elements[collect(domaincellset)]
    end
        
    return Grid(elements, nodes, facesets=facesets, cellsets=cellsets)
end

export gmsh
export tonodes, toelements, toboundary, tofacesets, tocellsets, saved_file_to_grid

end
