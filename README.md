# Ferrite.jl

<!---
[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg

[docs-dev-url]: http://koehlerson.github.io/gmsh.jl/dev/
-->

FerriteGmsh tries to simplify the conversion from a gmsh mesh to a Ferrite mesh.

## Installation

```
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

