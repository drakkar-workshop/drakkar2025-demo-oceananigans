# Drakkar 2025 Oceananigans / ClimaOcean

This repository contains the material used for the Oceananigans/ClimaOcean demonstration and hands-on-session @ Drakkar 2025 https://drakkar2025.sciencesconf.org

## Jupyter notebooks

The repository contains three different notebooks that will demostrate the capabilities of ClimaOcean and Oceananigans.

- _**cabbeling.ipnyb**_ shows how to use Oceananigans' `NonhydrostaticModel` to solve the Boussiesq Navier-Stokes equations in a two-dimensional x-z domain where density is calculated using the `TEOS10` equation of state.

- _**baroclinic_adjustment.ipynb**_ is adapted from the example in the Oceananigans repository (https://github.com/CliMA/Oceananigans.jl/blob/main/examples/baroclinic_adjustment.jl) and shows how to use the `HydrostaticFreeSurfaceModel` to solve the Hydrostatic approximation of the Boussinesq Navier Stokes equations with a free surface aka the **Primitive equations** in a simplified baroclinic instability test case.

- _**global_ocean_simulation.ipynb**_ shows how to set up a full-fledged global ocean simulation, that leverages [ClimaOcean.jl](https://github.com/CliMA/ClimaOcean.jl) and [OrthogonalSphericalShellGrids.jl](https://github.com/CliMA/OrthogonalSphericalShellGrids.jl) to run on a tripolar grid with a realistic bathymetry, initial conditions from ECCO climatology and (possibly) a realistic atmosphere from JRA55 reanalysis data. While the resolution is low to be able to allow running on a laptop, if a CUDA-enabled GPU is available, it is always possible leverage it to increase the resolution and see eddy activity develop.

---
_**NOTE:**_

The notebooks run comfortably on Windows / Linux / Mac machines. However, to exploit the real power of Oceananigans / ClimaOcean, it is recommended to obtain GPU access before the workshop to test the simulations on GPUs. At the moment only CUDA-enabled GPUs (NVIDIA) are compatible with Oceananigans (i.e., AMD and Intel GPUs will not work).

---

# Getting Started with Julia

Oceananigans.jl and ClimaOcean.jl are natively written in the [Julia](https://docs.julialang.org/en/v1/) language, so users should be familiar with the Julia syntax.

To download and install Julia for your system visit the [Julia download website](https://julialang.org/downloads/) (easy option) or follow the instructions to build Julia from source in the Julia [github repository](https://github.com/JuliaLang/julia) (a tad more difficult).

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
or by opening the downloaded julia binary.

The julia language comes with a native Package manager (Pkg) that allows downloading and installing registered Julia Packages directly through the REPL by doing

```julia
julia> using Pkg

julia> pkg"add MyPackage"
```


## Useful packages

In the notebooks, we use Julia packages other than the specific fluid dynamics / ocean modeling libraries:
- Oceananigans.jl
- ClimaOcean.jl
- OrthogonalSphericalShellGrids.jl

#### Jupyter notebooks: [IJulia.jl](https://github.com/JuliaLang/IJulia.jl)

**`IJulia.jl`** allows installing Julia kernels for usage in jupyter notebooks.
In addition, it is possible to customize the kernel with specific julia options, for example to install a multi-threaded kernel using 8 threads (and no optimization):

```julia
julia> using IJulia

julia> installkernel("Julia 8 threads", "--check-bounds=no", "-O0", "-t 8")
```

#### Plotting: [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie)

The native Julia plotting package for visualization