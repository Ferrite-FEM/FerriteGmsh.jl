# End-to-end tests on meshes that gmsh actually generates for the cell types added in
# https://github.com/Ferrite-FEM/FerriteGmsh.jl/issues/47: wedges from an extruded
# triangulation, and the hex/pyramid/tet mix gmsh produces where a recombined region
# meets a tetrahedral one.

"Total volume of `grid` and the smallest det(J) over all cells."
function grid_volume(grid)
    volume = 0.0
    mindetJ = Inf
    for cellid in 1:getncells(grid)
        cell = getcells(grid, cellid)
        ip = Ferrite.geometric_interpolation(typeof(cell))
        qr = QuadratureRule{Ferrite.getrefshape(ip)}(2)
        cv = CellValues(qr, ip, ip)
        reinit!(cv, getcoordinates(grid, cellid))
        for q in 1:getnquadpoints(cv)
            volume += getdetJdV(cv, q)
            mindetJ = min(mindetJ, getdetJdV(cv, q) / Ferrite.getweights(qr)[q])
        end
    end
    return volume, mindetJ
end

# Second order wedges ("Prism 18") need a Ferrite version that defines `QuadraticWedge`.
const WEDGE_ORDERS = isdefined(Ferrite, :QuadraticWedge) ?
    ((1, Ferrite.Wedge), (2, Ferrite.QuadraticWedge)) : ((1, Ferrite.Wedge),)

@testset "wedge mesh (order $order)" for (order, WedgeType) in WEDGE_ORDERS
    Gmsh.initialize()
    grid = try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("wedges")
        for (i, (x, y)) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)])
            gmsh.model.geo.addPoint(x, y, 0, 0.5, i)
        end
        for (i, (a, b)) in enumerate([(1, 2), (2, 3), (3, 4), (4, 1)])
            gmsh.model.geo.addLine(a, b, i)
        end
        gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
        gmsh.model.geo.addPlaneSurface([1], 1)
        # Extruding a triangulation with recombine=true keeps the layers as prisms
        # instead of subdividing them into tetrahedra.
        out = gmsh.model.geo.extrude([(2, 1)], 0, 0, 2, [3], Float64[], true)
        gmsh.model.geo.synchronize()
        # `extrude` returns the top surface first, then the volume, then the laterals.
        @assert out[1][1] == 2 && out[2][1] == 3
        topsurface = out[1][2]
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.setPhysicalName(2, 1, "bottom")
        gmsh.model.addPhysicalGroup(2, [topsurface], 2)
        gmsh.model.setPhysicalName(2, 2, "top")
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(order)
        togrid()
    finally
        Gmsh.finalize()
    end

    @test getncells(grid) > 0
    @test all(c -> c isa WedgeType, grid.cells)
    volume, mindetJ = grid_volume(grid)
    @test mindetJ > 0
    @test volume ≈ 2.0                              # 1 x 1 base extruded by 2

    # The triangular facet of a wedge is local facet 1, the top one facet 5.
    @test !isempty(getfacetset(grid, "bottom"))
    @test !isempty(getfacetset(grid, "top"))
    @test all(fi -> fi[2] == 1, getfacetset(grid, "bottom"))
    @test all(fi -> fi[2] == 5, getfacetset(grid, "top"))
    # Every triangle of the base surface must show up exactly once in each set.
    @test length(getfacetset(grid, "bottom")) == length(getfacetset(grid, "top"))
end

@testset "hex/pyramid/tet mesh" begin
    Gmsh.initialize()
    grid = try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("hybrid")
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
        gmsh.model.occ.addBox(0, 0, 1, 1, 1, 1, 2)
        gmsh.model.occ.fragment([(3, 1)], [(3, 2)])
        gmsh.model.occ.synchronize()
        # Structure and recombine the lower box into hexahedra; the upper box stays
        # tetrahedral, so gmsh fills the interface with pyramids.
        for (_, tag) in gmsh.model.getEntities(1)
            gmsh.model.mesh.setTransfiniteCurve(tag, 3)
        end
        for (_, tag) in gmsh.model.getBoundary([(3, 1)], false, false, false)
            gmsh.model.mesh.setTransfiniteSurface(abs(tag))
            gmsh.model.mesh.setRecombine(2, abs(tag))
        end
        gmsh.model.mesh.setTransfiniteVolume(1)
        gmsh.model.mesh.generate(3)
        togrid()
    finally
        Gmsh.finalize()
    end

    celltypes = Set(typeof.(grid.cells))
    @test Ferrite.Hexahedron in celltypes
    @test Ferrite.Pyramid in celltypes
    @test Ferrite.Tetrahedron in celltypes
    volume, mindetJ = grid_volume(grid)
    @test mindetJ > 0
    @test volume ≈ 2.0                              # two stacked unit boxes
end

@testset "boundary entity with mixed element types" begin
    # A single gmsh entity can hold elements of more than one type. Build a discrete
    # model where one surface entity carries both a triangle and a quadrilateral -- the
    # two kinds of facet a wedge has -- and check that neither is dropped.
    Gmsh.initialize()
    grid = try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("mixed_boundary")
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
                  (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0)]
        gmsh.model.addDiscreteEntity(3, 1)
        gmsh.model.mesh.addNodes(3, 1, collect(1:6), collect(Iterators.flatten(coords)))
        gmsh.model.mesh.addElements(3, 1, [6], [[1]], [collect(1:6)])   # one Prism 6
        gmsh.model.addDiscreteEntity(2, 1)
        # Triangle 3 = bottom face, Quadrilateral 4 = a lateral face of the same wedge.
        gmsh.model.mesh.addElements(2, 1, [2, 3], [[2], [3]], [[1, 3, 2], [1, 2, 5, 4]])
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.setPhysicalName(2, 1, "mixed")
        togrid()
    finally
        Gmsh.finalize()
    end

    @test getncells(grid) == 1
    # Both the triangular and the quadrilateral facet must be found.
    @test getfacetset(grid, "mixed") == Set([FacetIndex(1, 1), FacetIndex(1, 2)])
end
