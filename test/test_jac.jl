@testset "det J test" begin 
    ip = Serendipity{3,RefCube,2}()
    qr = QuadratureRule{3,RefCube}(2)
    cv = CellScalarValues(qr,ip)

    grid = togrid("quadhex_serendipity1.msh")
    dh = DofHandler(grid)
    push!(dh, :u, 2)
    close!(dh)

    for cell in CellIterator(dh)
        reinit!(cv,cell)
        @test all(cv.detJdV .≈ 0.125)
    end
end
