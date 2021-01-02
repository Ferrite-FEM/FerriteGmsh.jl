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

function getnodes()
    nodeid, nodes = JuAFEMGmsh.gmsh.model.mesh.getNodes()
    dim = Int64(JuAFEMGmsh.gmsh.model.getDimension()) #Int64 otherwise julia crashes
    return [Node(Vec{dim}(nodes[i:i+(dim-1)])) for i in 1:3:length(nodes)]
end

function getelements(dim)
    elementtypes, elementtags, nodetags = JuAFEMGmsh.gmsh.model.mesh.getElements(dim)
    @assert length(elementtypes)==1 "only one element type per mesh is supported"
    elementname, dim, order, numnodes, localnodecoord, numprimarynodes = JuAFEMGmsh.gmsh.model.mesh.getElementProperties(elementtypes[1]) 
    nodetags = convert(Array{Array{Int64,1},1}, nodetags)[1]
    juafemcell = gmshtojuafemcell[elementname]
    elements = [juafemcell(Tuple(nodetags[i:i+(numnodes-1)])) for i in 1:numnodes:length(nodetags)]
    return elements
end

function getboundary(dim)
    boundarydict = Dict{String, Vector}()
    boundaries = gmsh.model.getPhysicalGroups(dim) 
    for (tag,boundary) in enumerate(boundaries)
        physicaltag = boundary[2]
        name = gmsh.model.getPhysicalName(dim, physicaltag)
        boundarytypes, boundarytags, boundarynodetags = gmsh.model.mesh.getElements(dim, tag)
        elementname, dim, order, numnodes, localnodecoord, numprimarynodes = JuAFEMGmsh.gmsh.model.mesh.getElementProperties(boundarytypes[1]) 
        boundarynodetags = convert(Array{Array{Int64,1},1}, boundarynodetags)[1]
        boundaryconnectivity = [Tuple(boundarynodetags[i:i+(numnodes-1)]) for i in 1:numnodes:length(boundarynodetags)]
        boundarydict[name] = boundaryconnectivity
    end 
    return boundarydict
end

export gmsh
export getnodes, getelements, getboundary

end
