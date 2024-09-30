@testset "det J test" begin 
    if FerriteV1
        ip = Serendipity{RefHexahedron, 2}()
        qr = QuadratureRule{RefHexahedron}(2)
        cv = CellValues(qr, ip, ip)
    else
        ip = Serendipity{3,RefCube,2}()
        qr = QuadratureRule{3,RefCube}(2)
        cv = CellScalarValues(qr,ip)
    end

    grid = togrid("quadhex_serendipity1.msh")
    dh = DofHandler(grid)
    if FerriteV1
        add!(dh, :u, ip^2)
    else
        push!(dh, :u, 2)
    end
    close!(dh)

    for cell in CellIterator(dh)
        reinit!(cv, cell)
        @test all(cv.detJdV .≈ 0.125)
    end
end
