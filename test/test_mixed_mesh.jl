@testset "mixed mesh" begin
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("demo")

    lc = 0.2
    gmsh.model.geo.addPoint(-0.5, -1, 0, lc, 1)
    gmsh.model.geo.addPoint(0.5, -1, 0, lc, 2)
    gmsh.model.geo.addPoint(-0.5, 0, 0, lc, 3)
    gmsh.model.geo.addPoint(0.5, 0, 0, lc, 4)
    gmsh.model.geo.addPoint(-0.5, 1, 0, lc, 5)
    gmsh.model.geo.addPoint(0.5, 1, 0, lc, 6)

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 4, 2)
    gmsh.model.geo.addLine(4, 3, 3)
    gmsh.model.geo.addLine(1, 3, 4)
    gmsh.model.geo.addLine(3, 5, 5)
    gmsh.model.geo.addLine(5, 6, 6)
    gmsh.model.geo.addLine(4, 6, 7)

    gmsh.model.geo.addCurveLoop([1, 2, 3, -4], 1)
    gmsh.model.geo.addCurveLoop([-3, 7, -6, -5], 2)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.mesh.setTransfiniteCurve(1, 3)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, 3)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, 3)
    gmsh.model.geo.mesh.setTransfiniteCurve(4, 3)
    gmsh.model.geo.mesh.setTransfiniteCurve(5, 3)
    gmsh.model.geo.mesh.setTransfiniteCurve(6, 3)
    gmsh.model.geo.mesh.setTransfiniteCurve(7, 3)
    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setRecombine(2, 1)

    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.setPhysicalName(2, 1, "quad")

    gmsh.model.addPhysicalGroup(2, [2], 2)
    gmsh.model.setPhysicalName(2, 2, "triangle")

    gmsh.model.addPhysicalGroup(1, [6], 3)
    gmsh.model.setPhysicalName(1, 3, "top")

    gmsh.model.addPhysicalGroup(1, [1], 4)
    gmsh.model.setPhysicalName(1, 4, "bottom")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    nodes,_ = tonodes()
    elements, gmsh_eleidx = toelements(2)
    boundarydict = toboundary(1)
    facetsets = tofacetsets(boundarydict,elements)
    cellsets = tocellsets(2,gmsh_eleidx)
    if FerriteV1
        grid = Grid(elements,nodes,facetsets=facetsets,cellsets=cellsets)
    else
        grid = Grid(elements,nodes,facesets=facetsets,cellsets=cellsets)
    end
    @test getfacetset(grid, "bottom") == Set([FacetIndex(15,1),FacetIndex(17,1)])
    @test getfacetset(grid, "top") == Set([FacetIndex(3,3),FacetIndex(7,3)])
    @test grid.cellsets["triangle"] == Set([1,2,3,4,5,6,7,8,9,10,11,12,13,14])
    @test grid.cellsets["quad"] == Set([15,16,17,18])
    @test grid.nodes[12].x ≈ Vec{2}((0.0,1.0))
    @test grid.nodes[6].x ≈ Vec{2}((0.5,1.0))
    @test grid.nodes[5].x ≈ Vec{2}((-0.5,1.0))
    @test grid.nodes[7].x ≈ Vec{2}((0.0,-1.0))
    @test grid.nodes[2].x ≈ Vec{2}((0.5,-1.0))
    @test grid.nodes[1].x ≈ Vec{2}((-0.5,-1.0))
end
