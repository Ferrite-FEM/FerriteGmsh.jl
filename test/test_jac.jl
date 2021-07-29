using FerriteGmsh
using Test
using Ferrite
using Debugger

ip = Lagrange{3,RefCube,2}()
qr = QuadratureRule{3,RefCube}(2)
cv = CellScalarValues(qr,ip)

grid = saved_file_to_grid("test/quadhex_serendipity1.msh")
dh = DofHandler(grid)
push!(dh, :u, 2)
close!(dh)

for cell in CellIterator(dh)
    reinit!(cv,cell)
    @test all(cv.detJdV .≈ 0.125)
end
