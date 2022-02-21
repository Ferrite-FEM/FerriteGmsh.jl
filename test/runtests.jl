using FerriteGmsh
using Test
using Ferrite

include("test_jac.jl")
include("test_mixed_mesh.jl")
include("test_multiple_entities_group.jl")

@test_throws SystemError saved_file_to_grid("this-file-does-not-exist.msh")
