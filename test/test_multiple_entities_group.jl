@testset "multiple entities cellset" begin
    grid = saved_file_to_grid("matrix-inclusion.msh")
    @test grid.cellsets["inclusions"] == Set([5,6,182,164,153,186,185,4,13,168,183,
                                      158,176,177,179,12,173,171,11,162,7,172,156,
                                      2,10,167,157,169,180,160,9,161,159,152,170,163,
                                      8,184,1,154,175,166,14,3,178,155,181,15,174,165])
    @test grid.facesets["left"] == Set([Ferrite.FaceIndex((94, 3)),
                                      Ferrite.FaceIndex((113, 3)),
                                      Ferrite.FaceIndex((174, 3)),
                                      Ferrite.FaceIndex((87, 3)),
                                      Ferrite.FaceIndex((54, 3)),
                                      Ferrite.FaceIndex((134, 3)),
                                      Ferrite.FaceIndex((169, 2)),
                                      Ferrite.FaceIndex((175, 3)),
                                      Ferrite.FaceIndex((135, 3))])
    @test grid.facesets["right"] == Set([Ferrite.FaceIndex((92, 3)),
                                      Ferrite.FaceIndex((123, 3)),
                                      Ferrite.FaceIndex((136, 3)),
                                      Ferrite.FaceIndex((167, 3)),
                                      Ferrite.FaceIndex((164, 2)),
                                      Ferrite.FaceIndex((98, 3)),
                                      Ferrite.FaceIndex((62, 3)),
                                      Ferrite.FaceIndex((95, 3)),
                                      Ferrite.FaceIndex((165, 3))])
    @test grid.facesets["bottom"] == Set([Ferrite.FaceIndex((53, 3)),
                                        Ferrite.FaceIndex((89, 3)),
                                        Ferrite.FaceIndex((97, 3)),
                                        Ferrite.FaceIndex((156, 1)),
                                        Ferrite.FaceIndex((88, 3)),
                                        Ferrite.FaceIndex((159, 3)),
                                        Ferrite.FaceIndex((157, 3)),
                                        Ferrite.FaceIndex((112, 3)),
                                        Ferrite.FaceIndex((110, 3))])
    @test grid.facesets["top"] == Set([Ferrite.FaceIndex((124, 3)),
                                      Ferrite.FaceIndex((61, 3)),
                                      Ferrite.FaceIndex((184, 2)),
                                      Ferrite.FaceIndex((183, 3)),
                                      Ferrite.FaceIndex((114, 3)),
                                      Ferrite.FaceIndex((186, 3)),
                                      Ferrite.FaceIndex((91, 3)),
                                      Ferrite.FaceIndex((86, 3)),
                                      Ferrite.FaceIndex((93, 3))])
end
