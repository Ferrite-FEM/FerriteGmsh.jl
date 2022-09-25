using FerriteGmsh
using Test
using Ferrite

include("test_jac.jl")
include("test_mixed_mesh.jl")
include("test_multiple_entities_group.jl")
include("test_togrid.jl")

@test_throws SystemError togrid("this-file-does-not-exist.msh")
