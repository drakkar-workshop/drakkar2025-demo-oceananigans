# Drakkar 2025 Oceananigans / ClimaOcean 

This repository contains materials for the Oceananigans/ClimaOcean demonstration and hands-on session at [Drakkar 2025](https://drakkar2025.sciencesconf.org).

## Jupyter Notebooks

The repository includes three notebooks (and respective Julia scripts) showcasing the capabilities of ClimaOcean and Oceananigans:

1. **`cabbeling.ipynb`**  
   Demonstrates how to use Oceananigans' `NonhydrostaticModel` to solve the Boussinesq Navier-Stokes equations in a 2D x-z domain, with density computed using the `TEOS10` equation of state.

2. **`baroclinic_adjustment.ipynb`**  
   An adapted example from the [Oceananigans repository](https://github.com/CliMA/Oceananigans.jl/blob/main/examples/baroclinic_adjustment.jl). It illustrates the use of the `HydrostaticFreeSurfaceModel` to solve the hydrostatic approximation of the Boussinesq Navier-Stokes equations with a free surface (also known as the **Primitive Equations**) in a simplified baroclinic instability test case.

3. **`global_ocean_simulation.ipynb`**  
   Explains how to set up a global ocean simulation using [ClimaOcean.jl](https://github.com/CliMA/ClimaOcean.jl) and [OrthogonalSphericalShellGrids.jl](https://github.com/CliMA/OrthogonalSphericalShellGrids.jl). The simulation employs a tripolar grid with realistic bathymetry, initial conditions from ECCO climatology, and optional atmospheric data from the JRA55 reanalysis. The provided low-resolution setup runs on a laptop but can leverage CUDA-enabled GPUs for higher resolution.

---
 _**NOTE:**_

The notebooks run comfortably on Windows / Linux / Mac machines. However, to fully utilize the real power of Oceananigans / ClimaOcean, it is recommended to obtain GPU access before the workshop to test the simulations on GPUs. At the moment only CUDA-enabled GPUs (NVIDIA) are compatible with Oceananigans (i.e., AMD and Intel GPUs will not work).

---

# Getting Started with <img src="https://julialang.org/assets/infra/logo.svg" alt="drawing" width="80" />

Oceananigans.jl and ClimaOcean.jl are natively developed in the [Julia](https://docs.julialang.org/en/v1/) programming language. Therefore, users should be familiar with Julia syntax.

To install Julia, visit the [official download page](https://julialang.org/downloads/) (easy option) or follow the [source build instructions](https://github.com/JuliaLang/julia) on GitHub (slightly more complicated).

Once julia is downloaded and installed, it is possible to access the Julia interactive session or [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/#The-Julia-REPL) through a terminal
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
Alternatively, execute the downloaded Julia binary directly to access a REPL.

## Installing packages

Julia includes a native package manager (Pkg) for managing dependencies directly through the REPL. 

Start by installing the [IJulia.jl](https://github.com/JuliaLang/IJulia.jl) package, to use Julia with Jupyter notebooks:

```julia
julia> using Pkg

julia> pkg"add IJulia"
```

To install a custom kernel, such as one using 8 threads:

```julia
julia> using IJulia

julia> installkernel("Julia 8 threads", "--check-bounds=no", "-O0", "-t 8")
```

All other packages will be installed using Pkg directly in the notebooks / scripts. These are:

- [ClimaOcean.jl](https://github.com/CliMA/ClimaOcean.jl): the ocean model developed by the CliMA project 
- [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl): a ocean-flavored fluid dynamics library
- [SeawaterPolynomials.jl](https://github.com/CliMA/SeawaterPolynomials.jl): a package for seawater density computations
- [OrthogonalSphericalShellGrids.jl](https://github.com/CliMA/OrthogonalSphericalShellGrids.jl): a gridding package to generate grids
- [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie): the native Julia package for plotting and visualization

# Repository Contents
1. Notebooks
Stored in the `notebook` folder. To run a notebook, launch Jupyter Notebook from the command line:  
`% jupyter notebook`  
Then, select the desired notebook and choose the correct kernel on the top-right side of the screen.

2. Julia Scrips  
The `julia-scripts` folder contains equivalent Julia scripts for the experiments, for users that dislike notebooks and prefer running in an interactive Julia session. To execute the scripts in the Julia REPL:  
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

julia> include("julia-scripts/baroclinic_adjustment.jl")
```


### Visualization generated by the notebooks:

https://github.com/user-attachments/assets/c00d47ea-4b16-42c5-b559-4f6deb24a247


https://github.com/user-attachments/assets/1743450e-2d0a-4634-8495-fa0c72c181e4


https://github.com/user-attachments/assets/4b1d3f64-aa4f-448e-92f8-6e4c11ecba9f
