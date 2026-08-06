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

# do stuff to describe your gmsh model

dim = Int64(gmsh.model.getDimension())
facedim = dim - 1

# renumber the gmsh entities, such that the final used entities start by index 1
# this step is crucial, because Ferrite.jl uses the implicit entity index based on an array index
# and thus, need to start at 1
gmsh.model.mesh.renumberNodes()
gmsh.model.mesh.renumberElements()

# transfer the gmsh information
nodes = tonodes()
elements, gmsh_elementidx = toelements(dim)
# repair cells that Gmsh emitted with the opposite orientation; must happen before
# `tofacetsets` below, because flipping a cell renumbers its local facets
reorient!(elements, nodes)
cellsets = tocellsets(dim, gmsh_elementidx)

# "Domain" is the name of a PhysicalGroup and saves all cells that define the computational domain
domaincellset = cellsets["Domain"]
elements = elements[collect(domaincellset)]

boundarydict = toboundary(facedim)
facetsets = tofacetsets(boundarydict, elements)
gmsh.finalize()

using Ferrite
Grid(elements, nodes, facetsets=facetsets, cellsets=cellsets)
```

## Elements numbering & Supported elements

Ferrite might have a different element node-numbering scheme if compared to Gmsh. For correct portability of a mesh from Gmsh, it is important to transform Gmsh numbering in the one that Ferrite is expecting. Having the same numbering is important because the numbering provides information about which basis function is placed at each position, as well as the orientation of the elements. 

`FerriteGmsh` supports **every cell type that Ferrite defines**. Two `Dict`s in `src/FerriteGmsh.jl` describe the translation:

- `gmshtoferritecell` maps the Gmsh element name (as reported by `gmsh.model.mesh.getElementProperties`) to the Ferrite cell type,
- `gmshtoferriteperm` holds the node permutation for those element types whose numbering differs between the two codes, such that `ferrite_cell_nodes[i] == gmsh_element_nodes[perm[i]]`. Element types absent from it use the identity permutation.

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
where the numbers refers to the global nodes ids, which can be easily visualized in the following manner using the package [FerriteViz](https://github.com/Ferrite-FEM/FerriteViz.jl).


```julia
using Ferrite
using FerriteViz

grid =  generate_grid(QuadraticQuadrilateral,(1,1))

FerriteViz.wireframe(grid,markersize=14,strokewidth=20,textsize = 25, nodelabels=true,celllabels=true)
```
![Imgur](https://i.imgur.com/58OCFgo.png)

The Ferrite numbering `(1, 3, 9, 7, 2, 6, 8, 4, 5)` would have to match the numbering of Gmsh (see [gmsh docs](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering)).

In the particular case of the QuadraticQuadrilateral, the numbering used in `Ferrite` and the numbering used in Gmsh matches, so this element only needs an entry in `gmshtoferritecell` and no entry in `gmshtoferriteperm`.

If the numbering does not match, like for example in the `QuadraticTetrahedron`, an entry in `gmshtoferriteperm` specifies the correct numbering. Such an entry is found by matching the Gmsh local node coordinates -- the fifth return value of `gmsh.model.mesh.getElementProperties` -- against `Ferrite.reference_coordinates(Ferrite.geometric_interpolation(cell))`, up to the affine map between the two reference domains.

### Elements supported (summary):

All 14 cell types that `Ferrite` defines are supported:

| Gmsh element | Ferrite cell | reordered |
| --- | --- | --- |
| `Line 2` | `Line` | |
| `Line 3` | `QuadraticLine` | |
| `Triangle 3` | `Triangle` | |
| `Triangle 6` | `QuadraticTriangle` | |
| `Quadrilateral 4` | `Quadrilateral` | |
| `Quadrilateral 8` | `SerendipityQuadraticQuadrilateral` | |
| `Quadrilateral 9` | `QuadraticQuadrilateral` | |
| `Tetrahedron 4` | `Tetrahedron` | |
| `Tetrahedron 10` | `QuadraticTetrahedron` | ✓ |
| `Hexahedron 8` | `Hexahedron` | |
| `Hexahedron 20` | `SerendipityQuadraticHexahedron` | ✓ |
| `Hexahedron 27` | `QuadraticHexahedron` | ✓ |
| `Prism 6` | `Wedge` | |
| `Pyramid 5` | `Pyramid` | ✓ |

`Wedge`, `Pyramid` and `SerendipityQuadraticQuadrilateral` require Ferrite v1.

Gmsh element types that are not in this table -- the higher order variants such as `Triangle 10`, `Tetrahedron 20`, `Prism 18` or `Pyramid 14` -- have no counterpart in `Ferrite`, which only provides first and second order cells. Meshes containing them raise an error naming the offending element type; re-mesh with `gmsh.model.mesh.setOrder(1)` or `setOrder(2)`.

### Element orientation

Gmsh numbers the nodes of an element following the orientation of the entity it belongs to, so a surface whose normal points along `-z` produces clockwise elements while `Ferrite` requires a positive Jacobian determinant. `togrid` detects such cells and relabels them with `reorient!`, which keeps a cell on exactly the same physical region and only reverses its orientation. Correctly oriented cells are left untouched, so `ReverseMesh Surface{...};` in the `.geo` file is no longer needed.
