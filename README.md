# FerriteGmsh.jl

<!---
[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg

[docs-dev-url]: http://koehlerson.github.io/gmsh.jl/dev/
-->

FerriteGmsh tries to simplify the conversion from a gmsh mesh to a Ferrite mesh.

## Installation

In order to use this package, you need to install `gmsh.jl` first

```
]add https://github.com/koehlerson/gmsh.jl.git
]add https://github.com/koehlerson/FerriteGmsh.jl.git
```

## Example

![Imgur](https://i.imgur.com/eC2W4SZ.png)

The above example is taken from the 5th tutorial of [gmsh.jl](https://github.com/koehlerson/gmsh.jl).
In this tutorial `domain = "10"` corresponds to a specific `PhysicalGroup`, that gathers all cells of the computational domain. 

This package offers two workflows. The user can either load an already meshed file `.msh` or use it in an interactive way.

The first approach can be achieved by

```julia
using FerriteGmsh

saved_file_to_grid("t5.msh",domain="10")
```

while the latter is done by

```julia
using FerriteGmsh

gmsh.initialize()
dim = Int64(gmsh.model.getDimension())
facedim = dim - 1

# do stuff to prescribe your gmsh model

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

## Elements numberimg & Supported elements

Ferrite might have a different element node-numbering scheme if compared to Gmsh. For correct portability of a mesh from Gmsh, it is important to transform Gmsh numbering in the one that Ferrite is expecting. Having the same numbering is important because the numbering provides information about which basis function is placed at each position, as well as the orientation of the elements. 

By default `FerriteGmsh` supports all Ferrite elements in which the numbering is the same to the one used in Gmsh.

To check the numbering used in Ferrite, we could for example generate a grid with a single element, for a QuadraticQuadrilateral that would be:

```julia
using Ferrite
grid =  generate_grid(QuadraticQuadrilateral,(1,1))
```

Accessing the cells of that element:

```julia
julia> grid.cells
1-element Vector{QuadraticQuadrilateral}: 
QuadraticQuadrilateral((1, 3, 9, 7, 2, 6, 8, 4, 5))
```
where the numbers refers to the global nodes ids, which can be easily visualized in the following manner using the package [FerriteVis](https://github.com/koehlerson/FerriteVis.jl).


```julia
using Ferrite
using FerriteVis

grid =  generate_grid(QuadraticQuadrilateral,(1,1))

FerriteVis.wireframe(grid,markersize=14,strokewidth=20,textsize = 25, nodelabels=true,celllabels=true)
```
![Imgur](https://i.imgur.com/58OCFgo.png)

The Ferrite numbering `(1, 3, 9, 7, 2, 6, 8, 4, 5)` would have to match the numbering of Gmsh (see [gmsh docs](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering)).

In the particular case of the QuadraticQuadrilateral, the numbering used in `Ferrite` and the numbering used in Gmsh matches, then as it was anticipated this type of element would be supported through the default method for [translate_elements](https://github.com/koehlerson/FerriteGmsh.jl/blob/6682d9d4d95189f4799da19690b8ff0f18a9e177/src/FerriteGmsh.jl#L17-L19).


If the numbering does not match like for example in the `QuadraticTetrahedron` an specific method for the function [translate_elements](https://github.com/koehlerson/FerriteGmsh.jl/blob/6682d9d4d95189f4799da19690b8ff0f18a9e177/src/FerriteGmsh.jl#L21-L36) has to be created. In this method, the correct numbering is specified.

### Elements supported (summary):
With the elements supported by default and specific `translate_elements` methods (for QuadraticTetrahedron and 3D Serendipity), all the linear and quadratic elements available in `Ferrite` are already supported by `FerriteGmsh`.
