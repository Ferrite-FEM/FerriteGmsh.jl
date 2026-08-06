# Tests for the cell orientation repair of
# https://github.com/Ferrite-FEM/FerriteGmsh.jl/issues/15.
#
# Reuses `single_element_grid`, `cell_volume` and `REFERENCE_ELEMENTS` from
# test_cell_types.jl.

@testset "orientation tables" begin
    # Rederive `FerriteGmsh.cellorientation` from Ferrite's reference coordinates. The
    # orientation reversing symmetry is a reflection: swap the first two reference
    # coordinates, or negate the coordinate of a line.
    reflect(x::Vec{1}) = Vec{1}((-x[1],))
    reflect(x::Vec{2}) = Vec{2}((x[2], x[1]))
    reflect(x::Vec{3}) = Vec{3}((x[2], x[1], x[3]))

    "Determinant of the simplex spanned by the nodes `s` of `x`; positive if right handed."
    function simplexdet(x, s)
        e = ntuple(k -> x[s[k + 1]] - x[s[1]], length(s) - 1)
        return FerriteGmsh._det(e)
    end

    @test issetequal(keys(FerriteGmsh.cellorientation), values(FerriteGmsh.gmshtoferritecell))

    for (name, C) in sort!(collect(FerriteGmsh.gmshtoferritecell); by = first)
        @testset "$name" begin
            simplex, flip = FerriteGmsh.cellorientation[C]
            refcoords = Ferrite.reference_coordinates(Ferrite.geometric_interpolation(C))

            # `flip` must be exactly the relabelling induced by the reflection, which
            # determines it uniquely, and must be a permutation of all local nodes.
            @test length(flip) == length(refcoords)
            @test sort(collect(flip)) == collect(eachindex(refcoords))
            @test all(i -> refcoords[flip[i]] ≈ reflect(refcoords[i]), eachindex(refcoords))
            # Induced by a reflection, so applying it twice is the identity.
            @test ntuple(i -> flip[flip[i]], length(flip)) == ntuple(identity, length(flip))

            # `simplex` must span a positively oriented simplex of the reference element,
            # and the flip must turn it into a negatively oriented one.
            @test length(simplex) == length(first(refcoords)) + 1
            @test allunique(simplex)
            @test simplexdet(refcoords, simplex) > 0
            @test simplexdet(refcoords[collect(flip)], simplex) < 0
        end
    end
end

@testset "reorient! round trip" begin
    # Flipping a valid cell must produce the same physical region with the opposite
    # orientation, and `reorient!` must turn it back into the original cell.
    for (gmshtype, (coords, exactvolume)) in sort!(collect(REFERENCE_ELEMENTS); by = first)
        grid = single_element_grid(gmshtype, coords)
        cell = getcells(grid, 1)
        _, flip = FerriteGmsh.cellorientation[typeof(cell)]
        flipped = typeof(cell)(ntuple(i -> cell.nodes[flip[i]], length(cell.nodes)))

        @testset "$(nameof(typeof(cell)))" begin
            @test flipped != cell
            cells = [flipped]
            reorient!(cells, grid.nodes)
            @test only(cells) == cell
            # The repaired cell still covers the original element exactly.
            volume, mindetJ = cell_volume(Grid(cells, grid.nodes))
            @test mindetJ > 0
            @test volume ≈ exactvolume

            # An already correctly oriented cell must be left alone.
            untouched = [cell]
            reorient!(untouched, grid.nodes)
            @test only(untouched) == cell
        end
    end
end

"""
    mirrored_single_element(gmshtype, coords)

Build the single gmsh element `gmshtype` on a *mirrored* copy of `coords`. Reflecting the
geometry while keeping gmsh's local node order makes gmsh describe a negatively oriented
element, which is exactly what it does for a surface whose normal points the other way.

Returns `(nodes, rawcell, grid)`, where `rawcell` comes straight out of `toelements` with
no reorientation applied and `grid` is the result of the full `togrid` pipeline.
"""
function mirrored_single_element(gmshtype::Int, coords::Vector{NTuple{3, Float64}})
    Gmsh.initialize()
    return try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("mirrored_element")
        _, dim, _, numnodes, _, _ = gmsh.model.mesh.getElementProperties(gmshtype)
        dim = Int(dim)
        # Swap x and y; a line lives on the x axis alone and has to be flipped along it.
        mirror(c) = dim == 1 ? (-c[1], c[2], c[3]) : (c[2], c[1], c[3])
        gmsh.model.addDiscreteEntity(dim, 1)
        nodetags = collect(1:numnodes)
        gmsh.model.mesh.addNodes(dim, 1, nodetags, collect(Iterators.flatten(map(mirror, coords))))
        gmsh.model.mesh.addElements(dim, 1, [gmshtype], [[1]], [nodetags])
        gmsh.model.mesh.renumberNodes()
        gmsh.model.mesh.renumberElements()
        (tonodes(), only(first(toelements(dim))), togrid())
    finally
        Gmsh.finalize()
    end
end

@testset "clockwise mesh of every cell type" begin
    # For each supported element type, hand gmsh a negatively oriented element and check
    # that the full `togrid` pipeline turns it into a valid Ferrite cell covering the same
    # region. Mirroring only changes the sign of det(J), so the exact volume is unchanged.
    for (gmshtype, (coords, exactvolume)) in sort!(collect(REFERENCE_ELEMENTS); by = first)
        nodes, rawcell, grid = mirrored_single_element(gmshtype, coords)
        simplex, flip = FerriteGmsh.cellorientation[typeof(rawcell)]

        @testset "$(nameof(typeof(rawcell)))" begin
            # The element really is the pathological case: untouched by `reorient!` it
            # would hand Ferrite a negative Jacobian.
            @test FerriteGmsh._simplexdet(rawcell.nodes, nodes, simplex) < 0

            # `togrid` must have applied exactly the documented relabelling.
            @test getncells(grid) == 1
            cell = getcells(grid, 1)
            @test cell isa typeof(rawcell)
            @test cell == typeof(rawcell)(ntuple(i -> rawcell.nodes[flip[i]],
                                                 length(rawcell.nodes)))

            # ... and the result is a usable cell of the right size.
            volume, mindetJ = cell_volume(grid)
            @test mindetJ > 0
            @test volume ≈ exactvolume
        end
    end
end

@testset "clockwise surface mesh (issue #15)" begin
    # The geometry from the issue. Its curve loop runs such that the surface normal points
    # along -z, so gmsh emits every triangle clockwise.
    Gmsh.initialize()
    grid = try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("issue15")
        points = [(0, 0), (0, 10), (10, 10), (10, 20), (20, 20), (20, 5), (30, 5), (30, 0)]
        for (i, (x, y)) in enumerate(points)
            gmsh.model.geo.addPoint(x, y, 0, 3.0, i)
        end
        for (tag, (a, b)) in zip(1:8, [(1, 2), (2, 3), (3, 4), (4, 5),
                                       (5, 6), (6, 7), (7, 8), (8, 1)])
            gmsh.model.geo.addLine(a, b, tag)
        end
        gmsh.model.geo.addCurveLoop([2, 3, 4, 5, 6, 7, 8, 1], 1)
        gmsh.model.geo.addPlaneSurface([1], 1)
        gmsh.model.geo.synchronize()
        # Sanity check that this geometry really is the pathological one.
        @test gmsh.model.getNormal(1, [0.0, 0.0])[3] < 0
        gmsh.model.mesh.generate(2)
        togrid()
    finally
        Gmsh.finalize()
    end

    @test getncells(grid) > 0
    @test all(c -> c isa Ferrite.Triangle, grid.cells)

    # Every cell must now have a positive Jacobian, which is what Ferrite requires and
    # what used to fail. `reinit!` throws on a non-positive det(J), so this also covers
    # the case of `cell_volume` erroring out.
    volume = 0.0
    for cellid in 1:getncells(grid)
        cellvolume, mindetJ = cell_volume(Grid([getcells(grid, cellid)], grid.nodes))
        @test mindetJ > 0
        volume += cellvolume
    end
    # Area of the polygon by the shoelace formula.
    pts = [(0, 0), (0, 10), (10, 10), (10, 20), (20, 20), (20, 5), (30, 5), (30, 0)]
    exact = abs(sum(pts[i][1] * pts[mod1(i + 1, end)][2] -
                    pts[mod1(i + 1, end)][1] * pts[i][2] for i in eachindex(pts))) / 2
    @test volume ≈ exact
end

@testset "counter-clockwise surface mesh is untouched" begin
    # The same geometry with the curve loop traversed the other way already satisfies
    # Ferrite's convention, and must come out of `togrid` unchanged.
    Gmsh.initialize()
    nodes, elements = try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("ccw")
        points = [(0, 0), (0, 10), (10, 10), (10, 20), (20, 20), (20, 5), (30, 5), (30, 0)]
        for (i, (x, y)) in enumerate(points)
            gmsh.model.geo.addPoint(x, y, 0, 3.0, i)
        end
        for (tag, (a, b)) in zip(1:8, [(2, 1), (3, 2), (4, 3), (5, 4),
                                       (6, 5), (7, 6), (8, 7), (1, 8)])
            gmsh.model.geo.addLine(a, b, tag)
        end
        gmsh.model.geo.addCurveLoop([1, 8, 7, 6, 5, 4, 3, 2], 1)
        gmsh.model.geo.addPlaneSurface([1], 1)
        gmsh.model.geo.synchronize()
        @test gmsh.model.getNormal(1, [0.0, 0.0])[3] > 0
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.renumberNodes()
        gmsh.model.mesh.renumberElements()
        tonodes(), first(toelements(2))
    finally
        Gmsh.finalize()
    end

    before = copy(elements)
    reorient!(elements, nodes)
    @test elements == before
end
