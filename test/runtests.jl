using FerriteGmsh
using FerriteGmsh: Gmsh
using Test
using Ferrite

const FerriteV1 = isdefined(Ferrite, :FacetIndex)
if FerriteV1
    const getfacetset = Ferrite.getfacetset
    const FacetIndex = Ferrite.FacetIndex
else
    const getfacetset = Ferrite.getfaceset
    const FacetIndex = Ferrite.FaceIndex
end

include("test_jac.jl")
include("test_mixed_mesh.jl")
include("test_multiple_entities_group.jl")
include("test_togrid.jl")
include("test_saveall_flag.jl")
include("test_grid_sets_examples.jl")

@test_throws SystemError togrid("this-file-does-not-exist.msh")
