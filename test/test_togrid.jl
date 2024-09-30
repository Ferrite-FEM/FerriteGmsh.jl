using FerriteGmsh, Gmsh

@testset "togrid" begin
    Gmsh.initialize()
    # Points
    p1 = gmsh.model.geo.addPoint(0, 0, 0)
    p2 = gmsh.model.geo.addPoint(1, 0, 0)
    p3 = gmsh.model.geo.addPoint(1, 1, 0)
    p4 = gmsh.model.geo.addPoint(0, 1, 0)
    # Lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    # Surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surface = gmsh.model.geo.addPlaneSurface([loop])
    # Generate mesh
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate()
    # Fetch the grid
    grid = togrid()
    @test length(grid.cells) > 0 # smoke test...
    @test Bool(gmsh.isInitialized()) # Test that FerriteGmsh didn't finalize
    
    # gmsh.clear()
    gmsh.model.occ.addBox(0,0,0, 1,1,1, 1) # Smoke test to se eif it crashes
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate()
    grid = togrid()
    @test length(grid.cells) > 0 # smoke test...
    @test Bool(gmsh.isInitialized()) # Test that FerriteGmsh didn't finalize

    gmsh.finalize()
    @test !Bool(gmsh.isInitialized()) # Test that gmsh could shutdown
end
