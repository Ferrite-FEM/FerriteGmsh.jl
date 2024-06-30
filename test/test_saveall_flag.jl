@testset "saveall flag" begin
    # Initialize gmsh
    Gmsh.initialize()

    # Add the points
    h = 0.05
    o = gmsh.model.geo.add_point(0.0, 0.0, 0.0, h)
    p1 = gmsh.model.geo.add_point(0.5, 0.0, 0.0, h)
    p2 = gmsh.model.geo.add_point(1.0, 0.0, 0.0, h)
    p3 = gmsh.model.geo.add_point(0.0, 1.0, 0.0, h)
    p4 = gmsh.model.geo.add_point(0.0, 0.5, 0.0, h)

    # Add the lines
    l1 = gmsh.model.geo.add_line(p1, p2)
    l2 = gmsh.model.geo.add_circle_arc(p2, o, p3)
    l3 = gmsh.model.geo.add_line(p3, p4)
    l4 = gmsh.model.geo.add_circle_arc(p4, o, p1)

    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4])
    surf = gmsh.model.geo.add_plane_surface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "right")
    gmsh.model.add_physical_group(1, [l2], -1, "outer")
    gmsh.model.add_physical_group(1, [l3], -1, "left")
    gmsh.model.add_physical_group(1, [l4], -1, "inner")
    gmsh.model.add_physical_group(2, [surf])

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # Convert the mesh to Ferrite Grid
    grid_bad = FerriteGmsh.togrid()
    if FerriteV1
        close(VTKFile("grid_bad", grid_bad))
    else
        vtk_save(vtk_grid("grid_bad", grid_bad))
    end

    grid_good = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end
    @test gmsh.option.getNumber("Mesh.SaveAll") == 0
    @test Node((0.0,0.0)) ∉ grid_good.nodes
    nodeid = findall(x->x==Node((0.0,0.0)),grid_bad.nodes)
    @test !(any(findall(x->nodeid ∈ x.nodes, grid_bad.cells)))
end
