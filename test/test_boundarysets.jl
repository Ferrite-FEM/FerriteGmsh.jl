using FerriteGmsh
using Ferrite

grid = saved_file_to_grid("test/holed_plate.msh")

dh = DofHandler(grid)
push!(dh, :u, 1) 
close!(dh)

ch = ConstraintHandler(dh)
dbc_bottom = Dirichlet(:u, getfaceset(grid,"bottom"), (x,t) -> 0)
dbc_top = Dirichlet(:u, getfaceset(grid,"top"), (x,t) -> 1)
dbc_left = Dirichlet(:u, getfaceset(grid,"left"), (x,t) -> 2)
dbc_right = Dirichlet(:u, getfaceset(grid,"right"), (x,t) -> 3)
dbc_hole = Dirichlet(:u, getfaceset(grid,"hole"), (x,t) -> 4)
add!(ch,dbc_bottom)
add!(ch,dbc_top)
add!(ch,dbc_left)
add!(ch,dbc_right)
add!(ch,dbc_hole)
close!(ch)

vtk = vtk_grid("test_grid", grid)
vtk_point_data(vtk, ch)
vtk_save(vtk) 
