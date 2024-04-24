# FerriteGmsh.jl

<!---
[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg

[docs-dev-url]: http://ferrite-fem.github.io/FerriteGmsh.jl/dev/
-->

FerriteGmsh tries to simplify the conversion from a gmsh mesh to a Ferrite mesh.

## Installation

```
]add FerriteGmsh
```

## Example

![Imgur](https://i.imgur.com/qzQKx4x.png)

The above example is taken from the `test_mixed_grid.jl` file which can be found in the `test` folder.

This package offers two workflows. The user can either load an already defined geometry by a `.msh,.geo` file or use it in an interactive way.
The first approach can be achieved by

```julia
using FerriteGmsh

togrid("meshfile.msh")
togrid("meshfile.geo") #gets meshed automatically
```

while the latter is done by

```julia
using FerriteGmsh

gmsh.initialize()
dim = Int64(gmsh.model.getDimension())
facedim = dim - 1

# do stuff to describe your gmsh model

# renumber the gmsh entities, such that the final used entities start by index 1
# this step is crucial, because Ferrite.jl uses the implicit entity index based on an array index
# and thus, need to start at 1
gmsh.model.mesh.renumberNodes()
gmsh.model.mesh.renumberElements()

# transfer the gmsh information
nodes = tonodes()
elements, gmsh_elementidx = toelements(dim)
cellsets = tocellsets(dim, gmsh_elementidx)

# "Domain" is the name of a PhysicalGroup and saves all cells that define the computational domain
domaincellset = cellsets["Domain"]
elements = elements[collect(domaincellset)]

boundarydict = toboundary(facedim)
facesets = tofacesets(boundarydict, elements)
gmsh.finalize()

Grid(elements, nodes, facesets=facesets, cellsets=cellsets)
```

## Elements numbering & Supported elements

Ferrite might have a different element node-numbering scheme if compared to Gmsh. For correct portability of a mesh from Gmsh, it is important to transform Gmsh numbering in the one that Ferrite is expecting. Having the same numbering is important because the numbering provides information about which basis function is placed at each position, as well as the orientation of the elements. 

By default `FerriteGmsh` supports all Ferrite elements in which the numbering is the same to the one used in Gmsh.

To check the numbering used in Ferrite, we could for example generate a grid with a single element, for a QuadraticQuadrilateral that would be:

```julia
using Ferrite
grid = generate_grid(QuadraticQuadrilateral,(1,1))
```

Accessing the cells of that element:

```julia
julia> grid.cells
1-element Vector{QuadraticQuadrilateral}: 
QuadraticQuadrilateral((1, 3, 9, 7, 2, 6, 8, 4, 5))
```
where the numbers refers to the global nodes ids, which can be easily visualized in the following manner using the package [FerriteViz](https://github.com/koehlerson/FerriteViz.jl).


```julia
using Ferrite
using FerriteViz

grid =  generate_grid(QuadraticQuadrilateral,(1,1))

FerriteViz.wireframe(grid,markersize=14,strokewidth=20,textsize = 25, nodelabels=true,celllabels=true)
```
![Imgur](https://i.imgur.com/58OCFgo.png)

The Ferrite numbering `(1, 3, 9, 7, 2, 6, 8, 4, 5)` would have to match the numbering of Gmsh (see [gmsh docs](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering)).

In the particular case of the QuadraticQuadrilateral, the numbering used in `Ferrite` and the numbering used in Gmsh matches, then as it was anticipated this type of element would be supported through the default method for [translate_elements](https://github.com/koehlerson/FerriteGmsh.jl/blob/6682d9d4d95189f4799da19690b8ff0f18a9e177/src/FerriteGmsh.jl#L17-L19).


If the numbering does not match like for example in the `QuadraticTetrahedron` an specific method for the function [translate_elements](https://github.com/koehlerson/FerriteGmsh.jl/blob/6682d9d4d95189f4799da19690b8ff0f18a9e177/src/FerriteGmsh.jl#L21-L36) has to be created. In this method, the correct numbering is specified.

### Elements supported (summary):
With the elements supported by default and specific `translate_elements` methods (for QuadraticTetrahedron and 3D Serendipity), all the linear and quadratic elements available in `Ferrite` are already supported by `FerriteGmsh`.
