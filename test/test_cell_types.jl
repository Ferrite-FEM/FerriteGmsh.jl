# Pin down every entry of `FerriteGmsh.gmshtoferritecell` / `FerriteGmsh.gmshtoferriteperm`.
#
# For each supported gmsh element type we build a *single* element of a known, deliberately
# non-symmetric shape through the gmsh API, run it through `toelements`, and check that the
# resulting Ferrite cell describes the same geometry: the isoparametric map must have a
# positive Jacobian everywhere and the element must integrate to its exact volume. A wrong
# node permutation twists the element, which shows up as det(J) <= 0 or as a wrong volume.

"""
    single_element_grid(gmshtype, coords)

Build a mesh consisting of the single gmsh element of type `gmshtype` spanned by `coords`
(a vector of `NTuple{3,Float64}` in gmsh's local node order) and convert it with `togrid`.
"""
function single_element_grid(gmshtype::Int, coords::Vector{NTuple{3, Float64}})
    Gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("single_element")
        _, dim, _, numnodes, _, _ = gmsh.model.mesh.getElementProperties(gmshtype)
        @assert length(coords) == numnodes
        gmsh.model.addDiscreteEntity(dim, 1)
        nodetags = collect(1:numnodes)
        gmsh.model.mesh.addNodes(dim, 1, nodetags, collect(Iterators.flatten(coords)))
        gmsh.model.mesh.addElements(dim, 1, [gmshtype], [[1]], [nodetags])
        return togrid()
    finally
        Gmsh.finalize()
    end
end

"Integrate over cell 1 of `grid` using its geometric interpolation, returning (volume, min det(J))."
function cell_volume(grid)
    cell = getcells(grid, 1)
    ip = Ferrite.geometric_interpolation(typeof(cell))
    qr = QuadratureRule{Ferrite.getrefshape(ip)}(4)
    cv = CellValues(qr, ip, ip)
    reinit!(cv, getcoordinates(grid, 1))
    detJs = [getdetJdV(cv, q) / Ferrite.getweights(qr)[q] for q in 1:getnquadpoints(cv)]
    return sum(getdetJdV(cv, q) for q in 1:getnquadpoints(cv)), minimum(detJs)
end

# Reference geometries: gmsh element type => (nodes in gmsh local order, exact volume).
# The shapes are affine images of the reference element, so the exact volume is known and
# the quadrature is exact.
const REFERENCE_ELEMENTS = Dict{Int, Tuple{Vector{NTuple{3, Float64}}, Float64}}()

let
    # Affine maps chosen so that no two axes are scaled alike and the element is sheared,
    # which makes a wrong permutation impossible to hide.
    line(t) = (2.0 * t + 1.0, 0.0, 0.0)                       # [-1,1] -> [-1,3], length 4
    tri(u, v) = (2u + 0.5v, 3v, 0.0)                          # unit triangle, area 3
    quad(r, s) = (2r + 0.5s, 3s, 0.0)                         # [-1,1]^2, area 24
    tet(u, v, w) = (2u + 0.5v, 3v + 0.25w, 4w)                # unit tet, volume 24/6 = 4
    hex(r, s, t) = (2r + 0.5s, 3s + 0.25t, 4t)                # [-1,1]^3, volume 8*24 = 192
    prism(u, v, t) = (2u + 0.5v, 3v, 2t + 0.25v)              # tri x [-1,1], volume 3*2*2 = 12

    # Line 2 / Line 3
    REFERENCE_ELEMENTS[1] = ([line(-1.0), line(1.0)], 4.0)
    REFERENCE_ELEMENTS[8] = ([line(-1.0), line(1.0), line(0.0)], 4.0)
    # Triangle 3 / Triangle 6
    REFERENCE_ELEMENTS[2] = ([tri(0, 0), tri(1, 0), tri(0, 1)], 3.0)
    REFERENCE_ELEMENTS[9] = ([tri(0, 0), tri(1, 0), tri(0, 1),
                              tri(0.5, 0), tri(0.5, 0.5), tri(0, 0.5)], 3.0)
    # Quadrilateral 4 / 8 / 9
    let c = [quad(-1, -1), quad(1, -1), quad(1, 1), quad(-1, 1)],
        m = [quad(0, -1), quad(1, 0), quad(0, 1), quad(-1, 0)]
        REFERENCE_ELEMENTS[3] = (c, 24.0)
        REFERENCE_ELEMENTS[16] = (vcat(c, m), 24.0)
        REFERENCE_ELEMENTS[10] = (vcat(c, m, [quad(0, 0)]), 24.0)
    end
    # Tetrahedron 4 / 10
    REFERENCE_ELEMENTS[4] = ([tet(0, 0, 0), tet(1, 0, 0), tet(0, 1, 0), tet(0, 0, 1)], 4.0)
    REFERENCE_ELEMENTS[11] = ([tet(0, 0, 0), tet(1, 0, 0), tet(0, 1, 0), tet(0, 0, 1),
                               tet(0.5, 0, 0), tet(0.5, 0.5, 0), tet(0, 0.5, 0),
                               tet(0, 0, 0.5), tet(0, 0.5, 0.5), tet(0.5, 0, 0.5)], 4.0)
    # Hexahedron 8 / 20 / 27, in gmsh's local node order
    let c = [hex(-1, -1, -1), hex(1, -1, -1), hex(1, 1, -1), hex(-1, 1, -1),
             hex(-1, -1, 1), hex(1, -1, 1), hex(1, 1, 1), hex(-1, 1, 1)],
        e = [hex(0, -1, -1), hex(-1, 0, -1), hex(-1, -1, 0), hex(1, 0, -1),
             hex(1, -1, 0), hex(0, 1, -1), hex(1, 1, 0), hex(-1, 1, 0),
             hex(0, -1, 1), hex(-1, 0, 1), hex(1, 0, 1), hex(0, 1, 1)],
        f = [hex(0, 0, -1), hex(0, -1, 0), hex(-1, 0, 0), hex(1, 0, 0),
             hex(0, 1, 0), hex(0, 0, 1)]
        REFERENCE_ELEMENTS[5] = (c, 192.0)
        REFERENCE_ELEMENTS[17] = (vcat(c, e), 192.0)
        REFERENCE_ELEMENTS[12] = (vcat(c, e, f, [hex(0, 0, 0)]), 192.0)
    end
    # Prism 6
    REFERENCE_ELEMENTS[6] = ([prism(0, 0, -1), prism(1, 0, -1), prism(0, 1, -1),
                              prism(0, 0, 1), prism(1, 0, 1), prism(0, 1, 1)], 12.0)
    # Prism 18, in gmsh's local node order: vertices, the nine mid-edge nodes and the
    # centres of the three quadrilateral faces. Only for Ferrite versions defining
    # `QuadraticWedge`.
    if isdefined(Ferrite, :QuadraticWedge)
        REFERENCE_ELEMENTS[13] = ([prism(0, 0, -1), prism(1, 0, -1), prism(0, 1, -1),
                                   prism(0, 0, 1), prism(1, 0, 1), prism(0, 1, 1),
                                   prism(0.5, 0, -1), prism(0, 0.5, -1), prism(0, 0, 0),
                                   prism(0.5, 0.5, -1), prism(1, 0, 0), prism(0, 1, 0),
                                   prism(0.5, 0, 1), prism(0, 0.5, 1), prism(0.5, 0.5, 1),
                                   prism(0.5, 0, 0), prism(0, 0.5, 0), prism(0.5, 0.5, 0)], 12.0)
    end
    # Pyramid 5: 2x2 base in the z=0 plane (gmsh walks it counter-clockwise), apex sheared
    # off-centre so that a base rotation cannot pass by symmetry. V = 1/3 * 4 * 3 = 4.
    REFERENCE_ELEMENTS[7] = ([(-1.0, -1.0, 0.0), (1.0, -1.0, 0.0),
                              (1.0, 1.0, 0.0), (-1.0, 1.0, 0.0), (0.5, 0.25, 3.0)], 4.0)
end

@testset "cell types" begin
    # Every supported gmsh element type must be covered by a reference geometry, so that
    # adding an entry to `gmshtoferritecell` without a test fails here.
    supported = Dict{Int, String}()
    Gmsh.initialize()
    try
        for gmshtype in 1:19
            name, _, _, _, _, _ = gmsh.model.mesh.getElementProperties(gmshtype)
            haskey(FerriteGmsh.gmshtoferritecell, name) && (supported[gmshtype] = name)
        end
    finally
        Gmsh.finalize()
    end
    @test issetequal(keys(supported), keys(REFERENCE_ELEMENTS))

    for (gmshtype, name) in sort!(collect(supported))
        coords, exactvolume = REFERENCE_ELEMENTS[gmshtype]
        @testset "$name" begin
            grid = single_element_grid(gmshtype, coords)
            @test getncells(grid) == 1
            @test getcells(grid, 1) isa FerriteGmsh.gmshtoferritecell[name]
            volume, mindetJ = cell_volume(grid)
            @test mindetJ > 0
            @test volume ≈ exactvolume
        end
    end
end

@testset "unsupported cell type" begin
    # Triangle 10 (gmsh type 21) is a third order triangle; Ferrite has no such cell.
    Gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cubic_triangle")
        gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.setOrder(3)
        err = try togrid(); nothing catch e; e end
        @test err isa ErrorException
        @test occursin("unsupported gmsh element type \"Triangle 10\"", err.msg)
    finally
        Gmsh.finalize()
    end
end
