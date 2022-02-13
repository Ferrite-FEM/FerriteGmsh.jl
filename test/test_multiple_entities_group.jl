@testset "multiple entities cellset" begin
    grid = saved_file_to_grid("matrix-inclusion.msh")
    @test grid.cellsets["inclusions"] == Set([5,6,182,164,153,186,185,4,13,168,183,
                                      158,176,177,179,12,173,171,11,162,7,172,156,
                                      2,10,167,157,169,180,160,9,161,159,152,170,163,
                                      8,184,1,154,175,166,14,3,178,155,181,15,174,165])
end
