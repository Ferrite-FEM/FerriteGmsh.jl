module JuAFEMGmsh

using gmsh
using JuAFEM

const gmshtojuafemcell = Dict("Triangle 3" => Triangle)

function getnodes()
    nodeid, nodes = JuAFEMGmsh.gmsh.model.mesh.getNodes()
    dim = Int64(JuAFEMGmsh.gmsh.model.getDimension()) #Int64 otherwise julia crashes
    return [Node(Vec{dim}(nodes[i:i+(dim-1)])) for i in 1:3:length(nodes)]
end

function getelements(dim)
    elementtypes, elementtags, nodetags = JuAFEMGmsh.gmsh.model.mesh.getElements(dim)
    @assert length(elementtypes)==1 "only one element type per mesh is supported"
    elementname, dim, order, numnodes, localnodecoord, numprimarynodes = JuAFEMGmsh.gmsh.model.mesh.getElementProperties(2) 
    nodetags = convert(Array{Array{Int64,1},1}, nodetags)
    juafemcell = gmshtojuafemcell[elementname]
    elements = [juafemcell(Tuple(nodetags[1][i:i+(numnodes-1)])) for i in 1:numnodes:length(nodetags[1])]
    return elements
end

end
