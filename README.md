# JuAFEMGmsh.jl

<!---
[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg

[docs-dev-url]: http://koehlerson.github.io/gmsh.jl/dev/
-->

JuAFEMGmsh tries to simplify the conversion from a gmsh mesh to a JuAFEM mesh.

## Installation

```
]add https://github.com/koehlerson/JuAFEMGmsh.jl.git
```

## Example

![Imgur](https://imgur.com/eC2W4SZ)

The above example is taken from the 5th tutorial of [gmsh.jl](https://github.com/koehlerson/gmsh.jl).
In this tutorial `domain = "10"` corresponds to a specific `PhysicalGroup`, that gathers all cells of the computational domain. 

```julia
using JuAFEMGmsh

saved_file_to_grid("t5.msh",domain="10")
```

