using FerriteGmsh, Gmsh, LinearAlgebra

function generate_sample_data(dim)
    if dim == 3
        box = gmsh.model.occ.addBox(0,0,0,10,10,10)
        sphere = gmsh.model.occ.addSphere(5,5,5,2.5)
        res = gmsh.model.occ.fragment([(3,box)], [(3, sphere)])

        sphere = only(res[2][2])
        box = only(setdiff(Set(res[1]), Set([sphere])))

        gmsh.model.occ.synchronize()

        gmsh.model.addPhysicalGroup(3, [sphere[2]], 1, "inclusion")
        gmsh.model.addPhysicalGroup(3, [box[2]], 2, "matrix")

        # Add nodeset
        node_locations = getindex.(gmsh.model.getBoundingBox.(0, getindex.(gmsh.model.getEntities(0),2)),Ref(1:3))
        idx_origin = findfirst(==((0.0,0.0,0.0)),node_locations)
        idx_max = findfirst(==((10.0,10.0,10.0)),node_locations)
        gmsh.model.addPhysicalGroup(0, getindex.(gmsh.model.getEntities(0)[[idx_origin, idx_max]],2), 1, "extrema")

        # Add Faces
        for (i,(dim, face)) in enumerate(gmsh.model.getEntities(2))
            bb = gmsh.model.getBoundingBox(dim, face)
            close_to_min = isapprox.(bb,0.0,atol=1e-6)
            close_to_max = isapprox.(bb,10.0,atol=1e-6)
            if !any(close_to_min .|| close_to_max) ## inner surface
                gmsh.model.addPhysicalGroup(dim, [face], i, "inner surface")
                continue
            end
            zero_axis = (close_to_min[1:3] .&& close_to_min[4:end]) .|| (close_to_max[1:3] .&& close_to_max[4:end])
            if zero_axis[1]
                if close_to_min[1]
                    gmsh.model.addPhysicalGroup(dim, [face], i, "left")
                else
                    gmsh.model.addPhysicalGroup(dim, [face], i, "right")
                end
            elseif zero_axis[2]
                if close_to_min[2]
                    gmsh.model.addPhysicalGroup(dim, [face], i, "front")
                else
                    gmsh.model.addPhysicalGroup(dim, [face], i, "back")
                end
            else
                if close_to_min[3]
                    gmsh.model.addPhysicalGroup(dim, [face], i, "bottom")
                else
                    gmsh.model.addPhysicalGroup(dim, [face], i, "top")
                end
            end
        end

        gmsh.model.occ.synchronize()
    elseif dim ==2
        rectangle = gmsh.model.occ.addRectangle(0,0,0,10,10)
        circle = gmsh.model.occ.addDisk(5,5,0,2.5,2.5)
        res = gmsh.model.occ.fragment([(2,rectangle)], [(2, circle)])
        circle = only(res[2][2])
        rectangle = only(setdiff(Set(res[1]), Set([circle])))

        gmsh.model.occ.synchronize()

        gmsh.model.addPhysicalGroup(2, [circle[2]], 1, "inclusion")
        gmsh.model.addPhysicalGroup(2, [rectangle[2]], 2, "matrix")

        node_locations = getindex.(gmsh.model.getBoundingBox.(0, getindex.(gmsh.model.getEntities(0),2)),Ref(1:3))
        idx_origin = findfirst(==((0.0,0.0,0.0)),node_locations)
        idx_max = findfirst(==((10.0,10.0,0.0)),node_locations)
        gmsh.model.addPhysicalGroup(0, getindex.(gmsh.model.getEntities(0)[[idx_origin, idx_max]],2), 1, "extrema")

        for (i,(dim, face)) in enumerate(gmsh.model.getEntities(1))
            bb = gmsh.model.getBoundingBox(dim, face)[[1,2,4,5]]
            close_to_min = isapprox.(bb,0.0,atol=1e-6)
            close_to_max = isapprox.(bb,10.0,atol=1e-6)
            if !any(close_to_min .|| close_to_max) ## inner surface
                gmsh.model.addPhysicalGroup(dim, [face], i, "circular edge")
                continue
            end
            zero_axis = (close_to_min[1:2] .&& close_to_min[3:end]) .|| (close_to_max[1:2] .&& close_to_max[3:end])
            if zero_axis[1]
                if close_to_min[1]
                    gmsh.model.addPhysicalGroup(dim, [face], i, "left")
                else
                    gmsh.model.addPhysicalGroup(dim, [face], i, "right")
                end
            else
                if close_to_min[2]
                    gmsh.model.addPhysicalGroup(dim, [face], i, "bottom")
                else
                    gmsh.model.addPhysicalGroup(dim, [face], i, "top")
                end
            end
        end
    end
    mesh_size = 10.0
    gmsh.option.set_number("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.set_number("Mesh.MeshSizeMax", mesh_size)

    gmsh.model.mesh.generate(dim)
end

@testset "grid set examples" begin
    @testset "grid set example 2d" begin
        Gmsh.finalize()
        @assert Gmsh.initialize()
        dim = 2
        generate_sample_data(dim)

        gmsh.model.mesh.renumberNodes()
        gmsh.model.mesh.renumberElements()
        nodes, gmsh_nodeidx = tonodes()
        elements, gmsh_elementidx = toelements(dim)
        nodesets = FerriteGmsh.tonodesets(gmsh_nodeidx)
        cellsets = FerriteGmsh.tocellsets(dim, gmsh_elementidx)

        facets_boundary = toboundary(dim-1)
        facetsets = tofacetsets(facets_boundary, elements)

        grid = Grid(elements, nodes, facetsets=facetsets, cellsets=cellsets, nodesets=nodesets)

        # Check existance of sets
        @test "extrema" in keys(grid.nodesets)
        @test "top" in keys(grid.facetsets)
        @test "bottom" in keys(grid.facetsets)
        @test "left" in keys(grid.facetsets)
        @test "right" in keys(grid.facetsets)
        @test "circular edge" in keys(grid.facetsets)
        @test "inclusion" in keys(grid.cellsets)
        @test "matrix" in keys(grid.cellsets)

        # Get nodal coordinates of set extrema and compare to ferrite

        function get_facet_edgeset_coordinates(grid, setname)
            Set([get_node_coordinate(grid, node_id) for fidx in getfacetset(grid, setname) for node_id in Ferrite.edges(getcells(grid, fidx[1]))[fidx[2]]])
        end

        extrema_coordinates = get_node_coordinate.(getnodes(grid, "extrema"))
        @test any(isapprox(Ferrite.Vec{2}([0.0, 0.0])), extrema_coordinates)
        @test any(isapprox(Ferrite.Vec{2}([10.0, 10.0])), extrema_coordinates)
        @test length(extrema_coordinates) == 2

        edge_coords = get_facet_edgeset_coordinates(grid, "left")
        @test any(isapprox(Ferrite.Vec{2}([0.0, 0.0])), edge_coords)
        @test any(isapprox(Ferrite.Vec{2}([0.0, 10.0])), edge_coords)

        edge_coords = get_facet_edgeset_coordinates(grid, "right")
        @test any(isapprox(Ferrite.Vec{2}([10.0, 0.0])), edge_coords)
        @test any(isapprox(Ferrite.Vec{2}([10.0, 10.0])), edge_coords)

        edge_coords = get_facet_edgeset_coordinates(grid, "bottom")
        @test any(isapprox(Ferrite.Vec{2}([0.0, 0.0])), edge_coords)
        @test any(isapprox(Ferrite.Vec{2}([10.0, 0.0])), edge_coords)

        edge_coords = get_facet_edgeset_coordinates(grid, "top")
        @test any(isapprox(Ferrite.Vec{2}([0.0, 10.0])), edge_coords)
        @test any(isapprox(Ferrite.Vec{2}([10.0, 10.0])), edge_coords)

        # Check if nodes on inclusion boundary are in the correct distance to the center of the circle
        circular_edge_coords = get_facet_edgeset_coordinates(grid, "circular edge")
        @test all(map(collect(circular_edge_coords)) do x
            sum((x.-[5.0, 5.0]).^2) ≈ 2.5^2
        end)

        # Compute bounding box of inclusion
        inclusion_nodal_coordinates = [get_node_coordinate(grid, n) for c in getcells(grid, "inclusion") for n in c.nodes]
        x_ex, y_ex = extrema(getindex.(inclusion_nodal_coordinates,1)), extrema(getindex.(inclusion_nodal_coordinates,2))
        @test 2.5 <= x_ex[1] <= x_ex[2] <= 7.5
        @test 2.5 <= y_ex[1] <= y_ex[2] <= 7.5

        Gmsh.finalize()
    end
    @testset "grid set example 3d" begin
        Gmsh.finalize()
        @assert Gmsh.initialize()
        dim = 3
        generate_sample_data(dim)

        grid = togrid()

        # Check existance of sets
        @test "extrema" in keys(grid.nodesets)
        @test "top" in keys(grid.facetsets)
        @test "bottom" in keys(grid.facetsets)
        @test "left" in keys(grid.facetsets)
        @test "right" in keys(grid.facetsets)
        @test "front" in keys(grid.facetsets)
        @test "back" in keys(grid.facetsets)
        @test "inner surface" in keys(grid.facetsets)
        @test "inclusion" in keys(grid.cellsets)
        @test "matrix" in keys(grid.cellsets)

        # Get nodal coordinates of set extrema and compare to ferrite
        extrema_coordinates = get_node_coordinate.(getnodes(grid, "extrema"))
        @test any(isapprox(Ferrite.Vec{3}([0.0, 0.0, 0.0])), extrema_coordinates)
        @test any(isapprox(Ferrite.Vec{3}([10.0, 10.0, 10.0])), extrema_coordinates)
        @test length(extrema_coordinates) == 2

        # Check if all face normals point in the correct direction
        function compute_face_orientation(grid, setname)
            to_normal((n0, n1, n2)) = normalize((get_node_coordinate(grid, n1)-get_node_coordinate(grid, n0)) × (get_node_coordinate(grid, n2)-get_node_coordinate(grid, n0)))
            Set([to_normal(Ferrite.faces(getcells(grid, fidx[1]))[fidx[2]]) for fidx in getfacetset(grid, setname)])
        end

        face_normals = compute_face_orientation(grid, "left")
        @test all(x->x[1]≈-1, face_normals)
        face_normals = compute_face_orientation(grid, "right")
        @test all(x->x[1]≈1, face_normals)
        face_normals = compute_face_orientation(grid, "bottom")
        @test all(x->x[3]≈-1, face_normals)
        face_normals = compute_face_orientation(grid, "top")
        @test all(x->x[3]≈1, face_normals)
        face_normals = compute_face_orientation(grid, "back")
        @test all(x->x[2]≈1, face_normals)
        face_normals = compute_face_orientation(grid, "front")
        @test all(x->x[2]≈-1, face_normals)

        function get_facet_faceset_coordinates(grid, setname)
            Set([get_node_coordinate(grid, node_id) for fidx in getfacetset(grid, setname) for node_id in Ferrite.faces(getcells(grid, fidx[1]))[fidx[2]]])
        end

        # Check if nodes on inclusion boundary are in the correct distance to the center of the circle
        circular_edge_coords = get_facet_faceset_coordinates(grid, "inner surface")
        @test all(map(collect(circular_edge_coords)) do x
            sum((x.-[5.0, 5.0, 5.0]).^2) ≈ 2.5^2
        end)

        # Compute bounding box of inclusion
        inclusion_nodal_coordinates = [get_node_coordinate(grid, n) for c in getcells(grid, "inclusion") for n in c.nodes]
        x_ex, y_ex = extrema(getindex.(inclusion_nodal_coordinates,1)), extrema(getindex.(inclusion_nodal_coordinates,2))
        @test 2.5 <= x_ex[1] <= x_ex[2] <= 7.5
        @test 2.5 <= y_ex[1] <= y_ex[2] <= 7.5

        Gmsh.finalize()
    end

end