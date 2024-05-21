# ------------------------------------------------------------------------------
#
#  Gmsh Julia tutorial 16
#
#  Constructive Solid Geometry, OpenCASCADE geometry kernel
#
# ------------------------------------------------------------------------------

# Instead of constructing a model in a bottom-up fashion with Gmsh's built-in
# geometry kernel, starting with version 3 Gmsh allows you to directly use
# alternative geometry kernels. Here we will use the OpenCASCADE kernel.

Gmsh.initialize()

gmsh.model.add("t16")

# We first create two cubes:
gmsh.model.occ.addBox(0,0,0, 1,1,1, 1)
gmsh.model.occ.addBox(0,0,0, 0.5,0.5,0.5, 2)

# We apply a boolean difference to create the "cube minus one eigth" shape:
gmsh.model.occ.cut([(3,1)], [(3,2)], 3)

# Boolean operations with OpenCASCADE always create new entities. By default the
# extra arguments `removeObject' and `removeTool' in `cut()' are set to `True',
# which will delete the original entities.

# We then create the five spheres:
x = 0; y = 0.75; z = 0; r = 0.09

holes = []
for t in 1:5
    global x, z
    x += 0.166
    z += 0.166
    gmsh.model.occ.addSphere(x,y,z,r, 3 + t)
    t = (3, 3 + t)
    push!(holes, t)
end

# If we had wanted five empty holes we would have used `cut()' again. Here we
# want five spherical inclusions, whose mesh should be conformal with the mesh
# of the cube: we thus use `fragment()', which intersects all volumes in a
# conformal manner (without creating duplicate interfaces):
ov = gmsh.model.occ.fragment([(3,3)], holes)

gmsh.model.occ.synchronize()

lcar1 = .1
lcar2 = .0005
lcar3 = .055

# Assign a mesh size to all the points:
ov = gmsh.model.getEntities(0);
gmsh.model.mesh.setSize(ov, lcar1);

# Override this constraint on the points of the five spheres:
ov = gmsh.model.getBoundary(holes, false, false, true);
gmsh.model.mesh.setSize(ov, lcar3);

# Select the corner point by searching for it geometrically:
eps = 1e-3
ov = gmsh.model.getEntitiesInBoundingBox(0.5-eps, 0.5-eps, 0.5-eps,
                                         0.5+eps, 0.5+eps, 0.5+eps, 0)
gmsh.model.mesh.setSize(ov, lcar2)

gmsh.model.mesh.generate(3)

gmsh.write("t16.msh")
