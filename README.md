# Drakkar 2025 Oceananigans / ClimaOcean 

This repository contains the material for the Oceananigans/ClimaOcean demonstration and hands-on-session @ Drakkar 2025 https://drakkar2025.sciencesconf.org

## Jupyter notebooks

There are three different notebooks that demonstrate the capabilities of ClimaOcean and Oceananigans.

- _**cabbeling.ipynb**_ shows how to use Oceananigans' `NonhydrostaticModel` to solve the Boussinesq Navier-Stokes equations in a two-dimensional x-z domain where density is calculated using the [`TEOS10`](https://www.teos-10.org) equation of state.

- _**baroclinic_adjustment.ipynb**_ is an example adapted from the [Oceananigans documentation](https://clima.github.io/OceananigansDocumentation/stable/literated/baroclinic_adjustment/) that showcases how to use the `HydrostaticFreeSurfaceModel` to solve the hydrostatic approximation of the Boussinesq Navier Stokes equations with a free surface, also known as the **Primitive equations**, in a simplified baroclinic instability test case.

- _**global_ocean_simulation.ipynb**_ shows how to set up a fully-fledged global ocean simulation. It leverages [ClimaOcean.jl](https://github.com/CliMA/ClimaOcean.jl) and [OrthogonalSphericalShellGrids.jl](https://github.com/CliMA/OrthogonalSphericalShellGrids.jl) to run the simulation on a tripolar grid with a realistic bathymetry, initial conditions from ECCO climatology, and (possibly) a realistic atmosphere from JRA55 reanalysis data. While the resolution is low to allow running on a laptop, if a CUDA-enabled GPU is available, it is always possible to leverage it to increase the resolution and see eddy activity develop.

---
_**NOTE:**_

The notebooks run comfortably on Windows / Linux / Mac machines. However, to exploit the real power of Oceananigans / ClimaOcean, it is recommended to obtain GPU access before the workshop to test the simulations on GPUs. At the moment only CUDA-enabled GPUs (NVIDIA) are compatible with Oceananigans (i.e., AMD and Intel GPUs will not work).

---

# Getting Started with <img src="https://julialang.org/assets/infra/logo.svg" alt="drawing" width="80" />

Oceananigans.jl and ClimaOcean.jl are natively written in the [Julia](https://docs.julialang.org/en/v1/) language, so users should be familiar with the Julia syntax.

To download and install Julia for your system visit the [Julia download website](https://julialang.org/downloads/) (easy option) or follow the instructions to build Julia from source in the Julia [github repository](https://github.com/JuliaLang/julia) (a tad more difficult).

Once Julia is installed, we can access the Julia interactive session or [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julia-REPL) through a terminal
```julia
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.10.7 (2024-11-26)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```
or by opening the downloaded julia binary.

## Installing packages

Julia comes with a native Package manager (Pkg) that allows downloading and installing registered Julia Packages directly through the REPL.

The only package that needs to be installed is the [IJulia.jl](https://github.com/JuliaLang/IJulia.jl) package, that allows installing Julia kernels for usage in jupyter notebooks.

```julia
julia> using Pkg

julia> pkg"add IJulia"
```

Once **IJulia** is installed, it is possible to customize the kernel with specific julia options, for example to install a multi-threaded kernel using 8 threads (and no optimization):

```julia
julia> using IJulia

julia> installkernel("Julia 8 threads", "--check-bounds=no", "-O0", "-t 8")
```

All other packages will be installed using Pkg direcly in the notebooks. These are:

- [ClimaOcean.jl](https://github.com/CliMA/ClimaOcean.jl): the ocean model developed by the CliMA project 
- [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl): a ocean-flavored fluid dynamics library
- [SeawaterPolynomials.jl](https://github.com/CliMA/SeawaterPolynomials.jl): a package for seawater density computations
- [OrthogonalSphericalShellGrids.jl](https://github.com/CliMA/OrthogonalSphericalShellGrids.jl): a gridding package to generate grids
- [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie): the native Julia package for plotting and visualization

# Contents of the repository

The notebooks (in the `notebook` folder) can be executed by launching jupyter notebook via command line:  
`% jupyter notebook`  
After selecting the notebook, choose the correct kernel from the top right of the screen

The `julia-scripts` folder contain the same experiments written down as Julia scripts (for those who dislike notebooks). To run them in the REPL:
```julia
% julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.10.7 (2024-11-26)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> include("julia-scripts/baroclinic_adjustment.jl")
```


### Visualization generated by the notebooks:

https://github.com/user-attachments/assets/c00d47ea-4b16-42c5-b559-4f6deb24a247"


https://github.com/user-attachments/assets/1743450e-2d0a-4634-8495-fa0c72c181e4


https://github.com/user-attachments/assets/4b1d3f64-aa4f-448e-92f8-6e4c11ecba9f
