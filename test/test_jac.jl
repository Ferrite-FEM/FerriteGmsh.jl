@testset "det J test" begin 
    ip = Serendipity{RefHexahedron, 2}()
    qr = QuadratureRule{RefHexahedron}(2)
    cv = CellValues(qr, ip, ip)

    grid = togrid("quadhex_serendipity1.msh")
    dh = DofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    for cell in CellIterator(dh)
        reinit!(cv, cell)
        @test all(cv.detJdV .≈ 0.125)
    end
end
