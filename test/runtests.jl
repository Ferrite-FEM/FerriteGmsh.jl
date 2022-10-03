using FerriteGmsh
using FerriteGmsh: Gmsh
using Test
using Ferrite

include("test_jac.jl")
include("test_mixed_mesh.jl")
include("test_multiple_entities_group.jl")
include("test_togrid.jl")
include("test_saveall_flag.jl")

@test_throws SystemError togrid("this-file-does-not-exist.msh")
