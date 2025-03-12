
@testset "fix ordering" begin
    function randomize_cell_order(cell::CT) where {CT <: Ferrite.AbstractCell}
        nv = Ferrite.nvertices(cell)
        shift = rand(0:(nv - 1))
        flip = rand(Bool)
        if flip
            vertex_indices = ntuple(i -> 1 + nv - i, nv)
        else
            vertex_indices = ntuple(identity, nv)
        end
        vertex_indices = map(j -> mod1(j + shift, nv), vertex_indices)
        vertex_nodes = map(i -> cell.nodes[i], vertex_indices)
        if nv == length(cell.nodes)
            return CT(vertex_nodes)
        end
        edge_indices = ntuple(i -> findfirst(Ferrite.edges(cell)) do e 
            n1 = cell.nodes[vertex_indices[mod1(i, nv)]]
            n2 = cell.nodes[vertex_indices[mod1(i + 1, nv)]]
            e == (n1, n2) || e == (n2, n1)
        end, nv)
        edge_nodes = map(i -> cell.nodes[i + nv], edge_indices)
        if 2nv == length(cell.nodes)
            return CT((vertex_nodes..., edge_nodes...))
        else
            return CT((vertex_nodes..., edge_nodes..., cell.nodes[end]))
        end
    end

    function get_nodeset_from_set(grid, set::AbstractSet, entity_fun::Function)
        nodes = sizehint!(Set{Int}(), length(set))
        for (cellnr, entitynr) in set
            for node in entity_fun(getcells(grid, cellnr))[entitynr]
                push!(nodes, node)
            end
        end
        return nodes
    end
    function get_nodeset_from_set(grid::Grid, set::AbstractSet{FacetIndex})
        return get_nodeset_from_set(grid, set, Ferrite.facets)
    end
    function get_nodeset_from_set(grid::Grid, set::AbstractSet{VertexIndex})
        return get_nodeset_from_set(grid, set, Ferrite.vertices)
    end
    function export_nodes(vtk::VTKGridFile, grid, nodes, name)
        data = zeros(getnnodes(grid))
        for n in nodes
            data[n] = 1
        end
        write_node_data(vtk, data, name)
    end

    @testset "reverse_vertex_facet_mapping" begin
        t = Triangle((2, 4, 3))
        qt = QuadraticTriangle((2, 5, 3, 1, 10, 4))
        q = Quadrilateral((1, 3, 2, 5))
        qq = QuadraticQuadrilateral((4, 2, 1, 10, 5, 3, 11, 9, 8))
        sqq = SerendipityQuadraticQuadrilateral((1, 2, 3, 5, 4, 9, 8, 7))
        for c in (t, qt, q, qq, sqq)
            vertex_mapping = FerriteGmsh.reverse_vertex_mapping(c)
            facet_mapping = FerriteGmsh.reverse_facet_mapping(c)
            rc = FerriteGmsh.reverse_cell(c)
            @testset "vertices" begin
                for i in 1:Ferrite.nvertices(c)
                    @test Ferrite.vertices(c)[i] == Ferrite.vertices(rc)[vertex_mapping[i]]
                end
            end
            @testset "facets" begin
                for i in 1:nfacets(c)
                    @test Ferrite.facets(c)[i] == reverse(Ferrite.facets(rc)[facet_mapping[i]])
                end
            end
        end
    end

    function test_positive_detJ(grid::Grid{2, CT}) where {CT}
        @assert isconcretetype(CT)
        ip = geometric_interpolation(CT)
        RefShape = getrefshape(CT)
        qr = QuadratureRule{RefShape}(3)
        cv = CellValues(qr, ip, ip)
        for cell in CellIterator(grid)
            @test_nowarn reinit!(cv, cell)
        end
    end

    @testset "fix_node_order!(::Grid{2, $(CT)})" for CT in (Triangle, QuadraticTriangle, Quadrilateral, QuadraticQuadrilateral)
        grid = generate_grid(CT, (3, 3))
        map!(randomize_cell_order, grid.cells, grid.cells)
        empty!(grid.facetsets)

        for (key, component, value) in (("left", 1, -1), ("right", 1, 1), ("bottom", 2, -1), ("top", 2, 1))
            addfacetset!(grid, key, x -> x[component] ≈ value)
            addvertexset!(grid, key, x -> x[component] ≈ value)
        end

        grid_original = deepcopy(grid)
        FerriteGmsh.fix_node_order!(grid)
        test_positive_detJ(grid)
        for key in keys(grid.facetsets)
            nodes = get_nodeset_from_set(grid, getfacetset(grid, key))
            nodes_original = get_nodeset_from_set(grid_original, getfacetset(grid_original, key))
            @test nodes == nodes_original
        end
        for key in keys(grid.vertexsets)
            nodes = get_nodeset_from_set(grid, getvertexset(grid, key))
            nodes_original = get_nodeset_from_set(grid_original, getvertexset(grid_original, key))
            @test nodes == nodes_original
        end
    end
end
