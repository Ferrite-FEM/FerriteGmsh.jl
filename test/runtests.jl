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
# Wedge, Pyramid and SerendipityQuadraticQuadrilateral only exist on Ferrite v1, and so
# does the `geometric_interpolation`/`CellValues` API these tests are written against.
FerriteV1 && include("test_cell_types.jl")
FerriteV1 && include("test_hybrid_mesh.jl")
include("test_mixed_mesh.jl")
include("test_multiple_entities_group.jl")
include("test_togrid.jl")
include("test_saveall_flag.jl")

@test_throws SystemError togrid("this-file-does-not-exist.msh")
