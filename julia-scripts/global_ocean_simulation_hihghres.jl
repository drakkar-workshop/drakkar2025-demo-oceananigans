# # Realistic Ocean Simulations in pure Julia
#
# In this tutorial we will build a global ocean simulation using the ClimaOcean and Oceananigans Julia packages.
# The simulation will be at a nominal four-degree with parameterizations for mesoscale eddies and vertical mixing.
# We will set up the grid, download and interpolate the bathymetry, 
# download forcing files, and finally set up and run the simulation

# ### ***Bonus***
# At the end of the tutorial we will change the resolution to allow the spontaneous generation of eddies and 
# remove the eddy parameterization to see some beautiful ocean turbulence develop! (make sure to have GPU access
# for this step!)

# ## Required software
#
# The tutorial is quite computationally expensive, therefore, if you have 
# access, it is recommended to run the tutorial on one GPU. 
# However, for the purpose of understanding how the library works, 
# a 4 degree global ocean model runs comfortably on a laptop.

# ## Packages:

# Let's start by importing the necessary packages, these are:
# - ClimaOcean: the ocean model
# - Oceananigans: the fluid dynamics library doing the heavy lifting
# - OrthogonalSphericalShellGrids: contains the constructor for the Tripolar grid we will use in the tutorial
# - Printf: always useful for spitting out output
# - CairoMakie: visualization package to visualize the results

using Pkg
Pkg.activate("./")
pkg"add Oceananigans#ss/for-drakkar" 
pkg"add ClimaOcean@0.3.3"
pkg"add OrthogonalSphericalShellGrids"
pkg"add CairoMakie"
pkg"add CFTime"

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids
using OrthogonalSphericalShellGrids
using CFTime
using Dates
using Printf
using CairoMakie

# # Building a Global Ocean domain

# We will start building a global ocean in steps: 
# - (1) specifying an **Architecture**
# - (2) choosing a **Vertical coordinate** 
# - (2) Building a **Grid** 
# - (3) Downloading an interpolating a **Bathymetry**

# ### Architectures

# using an architecture is easy...
# it is possible to choose between:
# - (1) CPU 
# - (2) GPU 
# - (3) Distributed 

# In this case, we select `CPU()`, which is always the right choice to start prototyping. Building the script on CPUs is better to catch potential bugs in the script. 
# Once we know that everything works correctly, we can just change the following line to `arch = GPU()`, and voilá, the simulation runs on GPUs

arch = GPU()

# ### Vertical coordinates

# Oceananigans currently supports only $z$ and $z^\star$ coordinates.
# ClimaOcean provides a simple utility to build a simple exponential vertical coordinate

depth = 5000meters
Nz    = 50
h     = 10 

r_faces = ClimaOcean.exponential_z_faces(; Nz, h, depth)

# To use z-star coordinates we need to use a `z_faces = MutableVerticalDiscetization(r_faces)`, as opposed to a to 
# a `StaticVerticalDiscretization`, to set up the data structures required for a free-surface 
# following vertical coordinate. 

z_faces = MutableVerticalDiscretization(r_faces)

# ## Building a grid

# ClimaOcean allows building ocean simulations on three different grid types:
# - `RectilinearGrid`s which represent a _box_ or a Cartesian domain
# - `LatitudeLongitudeGrid`s, which discretizes the sphere along latitude and longitude lines
# - `OrthogonalSphericalShellGrid`s which discretize the sphere with two-dimensional coordinates that do not need to follow latitude and longitude lines. 
#    The only constraint is that the grid must be locally orthogonal in the horizontal direction.
#
# `LatitudeLongitudeGrid`s are the easiest grids to work with since coordinates are one-dimensional and globally orthogonal.
# (i.e. latitude depends only on the `j`-index and longitude depends only on the `i`-index)
# However, `LatitudeLongitudeGrid`s have the problem of the zonal spacing approaching zero as we move to the poles. 
#
# For this reason, to represent a global ocean we use a specific type of `OrthogonalSphericalShellGrid`, 
# called `TripolarGrid` that discretizes the sphere as a set of perpendicular ellipses and hyperbolae.
#
# Let's build a coarse `TripolarGrid` (about 4 degree resolution). 
# We pass to the grid, the architecture, the floating point precision, the size of the grid, and the vertical coordinate.

Nx = 1440 # longitudinal direction -> 1440 points is about 0.25ᵒ resolution
Ny = 650  # meridional direction -> same thing, 48 points is about 0.25ᵒ resolution
Nz = length(r_faces) - 1
grid = TripolarGrid(arch, Float64; size=(Nx, Ny, Nz), z=z_faces, halo=(7, 7, 4))

# ## Adding a bathymetry to the grid
#
# ClimaOcean provides a nifty utility to regrid the bathymetry over the grid, the `regrid_bathymetry` function.
# By default ClimaOcean downloads the ETOPO_2022 bathymetry at 1/60ᵒ resolution (459 MB) from the NOAA servers.
# However, since the servers are quite busy, I have uploaded a lower resolution version file to dropbox.
# !!! NOTE: This will download the ETOPO_2022 bathymetry, so make sure that you have an internet connection

url = "https://www.dropbox.com/scl/fi/zy1cu64ybn93l67rjgiq0/Downsampled_ETOPO_2022.nc?rlkey=5upqvoxrnljj205amqf663vcw&st=ou8b32tt&dl=0"
filename = isfile("Downsampled_ETOPO_2022.nc") ? "Downsampled_ETOPO_2022.nc" : download(url, "Downsampled_ETOPO_2022.nc")
bottom_height = regrid_bathymetry(grid; minimum_depth=10, major_basins=1, filename, dir="./")

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

# # Configuring an Ocean model
#
# ### Numerical details
#
# Oceananigans allows several numerical options. 
# We use a WENO schemes for the advection of momentum and 
# a centered scheme for tracer advection, to avoid implicit diapycnal diffusion of tracers.
# Stability in the momentum field is ensured by the WENO method. For the tracer field, since the centered
# scheme is dispersive, we need to add some explicit diffusion to avoid numerical instabilities.

momentum_advection = WENOVectorInvariant(order=5) 
tracer_advection   = WENO(order=5)

free_surface = SplitExplicitFreeSurface(grid; substeps=75) 

# ### Physical parameterizations
#
# We add a GM parameterization for mesoscale eddies and a CATKE vertical mixing scheme.
# All the closures require passing also the desired floating point precision of the model

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength

vertical_mixing = CATKEVerticalDiffusivity(mixing_length=CATKEMixingLength(Cᵇ = 0.001)) 
closure = vertical_mixing

# ### Building the ocean simulation
# 
# ClimaOcean provides a utility to build an ocean simulation with all the necessary components. 
# The function `ocean_simulation` returns a `Simulation` object of a `HydrostaticFreeSurfaceModel` that has 
# all the necessary components (BC, drag, etc) to run a global ocean simulation.

ocean = ocean_simulation(grid; 
                         momentum_advection, 
                         tracer_advection, 
                         closure, 
                         free_surface)

# ### Initialize our Ocean
#
# We use ECCO climatology to initialize the temperature and salinity fields. 
# We can use the metadata we defined earlier to set the initial conditions. 
 
temperature = ECCOMetadata(:temperature;  dir="./", version=ClimaOcean.ECCO.ECCO2Daily())
salinity    = ECCOMetadata(:salinity;     dir="./", version=ClimaOcean.ECCO.ECCO2Daily())

set!(ocean.model, T=temperature, S=salinity)

# # Adding an atmosphere
#
# ClimaOcean is a prototype for a coupled earth system model. 
# It couples an atmosphere to an ocean and computes air-sea fluxes using bulk formulae.
# At the moment, ClimaOcean provides a utility to download the JRA55 atmospheric reanalysis
# and use it as a prescribed atmosphere.
#
# !!! NOTE: This will download the JRA55 atmospheric reanalysis, so make sure that you have an internet connection (and enough disk space)
#
# We use an idealized atmosphere for this tutorial to avoid downloading the JRA55 data (~15GB).

atmosphere  = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(10))

# If we had a realistic atmosphere we would add radiative properties, however, since we do not have 
# downwelling radiation, we set it to nothing

radiation = Radiation(ocean_albedo = LatitudeDependentAlbedo())

# ### Coupling the atmosphere to the ocean
#
# The `OceanSeaIceModel` is an `AbstractModel` defined in ClimaOcean that couples an ocean to an atmosphere and a sea ice component.
# For the moment, the sea-ice component is not implemented, so we will only couple the ocean to the atmosphere.
# Instead of the sea ice model, we limit the temperature of the ocean to the freezing temperature.

similarity_theory = SimilarityTheoryTurbulentFluxes(grid)
sea_ice = ClimaOcean.FreezingLimitedOceanTemperature()
earth_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation, similarity_theory)

# ### Building the simulation
#
# We build the simulation out of the `earth_model` as we would do for any other Oceananigans model.
# We start with a smallish time-step (5 minutes) and run only for 2 days to dissipate initialization shocks.

earth = Simulation(earth_model; Δt=2minutes, stop_time=2days)

# ### Adding some diagnostics
#
# We add a callback to save surface fields as well as surface fluxes, every 6 hours

u, v, _ = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
s = sqrt(u^2 + v^2)

η = ocean.model.free_surface.η 

earth.output_writers[:surface_tracers] = JLD2OutputWriter(ocean.model, (; T, S, s),
                                                          schedule = TimeInterval(12hours),
                                                          indices = (:, :, grid.Nz),
                                                          overwrite_existing = true,
                                                          filename = "surface_fields.jld2")


earth.output_writers[:free_surface] = JLD2OutputWriter(ocean.model, (; η),
                                                       schedule = TimeInterval(12hours),
                                                       overwrite_existing = true,
                                                       filename = "free_surface.jld2")

Q  = earth.model.fluxes.total.ocean.heat
τx = earth.model.fluxes.total.ocean.momentum.u
τy = earth.model.fluxes.total.ocean.momentum.v
PE = earth.model.fluxes.total.ocean.tracers.S

earth.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, (; Q, τx, τy, PE),
                                                 schedule = TimeInterval(12hours),
                                                 overwrite_existing = true,
                                                 filename = "surface_fluxes.jld2")

# Also, we add a callback to print a message about how the simulation is going

wall_time = [time_ns()]

function progress(earth)
    clock   = earth.model.clock

    maxu = maximum(abs, u)
    maxv = maximum(abs, v)
    maxT = maximum(T)
    minS = minimum(S)
    
    @info @sprintf("Iteration: %d, time: %s, wall_time: %s, max(|u|, |v|): %.2e %.2e max(T): %.2e, min(S): %.2e\n",
                   clock.iteration, prettytime(clock.time), prettytime(1e-9 * (time_ns() - wall_time[1])), maxu, maxv, maxT, minS)

    wall_time[1] = time_ns()
end

add_callback!(earth, progress, IterationInterval(10))

# ### Running the simulation
#
# quite simply

run!(earth)

# Finished the two days, we increase timestep size and stop time and restart

earth.stop_time = 60days
earth.Δt = 10minutes

run!(earth)

# ## Visualizing the results
#
# We can visualize the results using CairoMakie. We record a video of surface variables and fluxes.
# To load the data we can use Oceananigans' `FieldTimeSeries` object.

using JLD2
using Oceananigans
using Oceananigans.Grids: halo_size
using CairoMakie 
using Statistics: mean

file  = jldopen("free_surface.jld2")
iters = keys(file["timeseries/t"]) 

Hx, Hy, _ = halo_size(η.grid)
T  = FieldTimeSeries("surface_fields.jld2", "T")
S  = FieldTimeSeries("surface_fields.jld2", "S")
s  = FieldTimeSeries("surface_fields.jld2", "s")

n  = Observable(1)
Tn = @lift(interior(T[$n], :, :, 1))
Sn = @lift(interior(S[$n], :, :, 1))
sn = @lift(interior(s[$n], :, :, 1))
ηn = @lift(file["timeseries/η/" * iters[$n]][Hx+1:end-Hx, Hy+1:end-Hy, 1])

fig = Figure(size = (1800, 800))
axT = Axis(fig[1, 1], title="Surface temperature ᵒC")
axS = Axis(fig[1, 2], title="Surface salinity psu")
axs = Axis(fig[2, 1], title="Surface speed ms⁻¹")
axη = Axis(fig[2, 2], title="Sea surface height m")

λ, φ, z = nodes(T[1])

hmT = heatmap!(axT, Tn, colormap=:magma,  colorrange=(-1, 30))
hmS = heatmap!(axS, Sn, colormap=:haline, colorrange=(25, 40))
hms = heatmap!(axs, sn, colormap=:deep,   colorrange=( 0, 0.8))
hmη = heatmap!(axη, ηn, colormap=:bwr,    colorrange=(-1, 1))

CairoMakie.record(fig, "surface_fields.mp4", 1:length(T.times); framerate=5) do i 
    @info "doing $i of $(length(T.times))"
    n[] = i
end

# let's also visualize the surface fluxes that force the model

Q  = FieldTimeSeries("surface_fluxes.jld2", "Q")
τx = FieldTimeSeries("surface_fluxes.jld2", "τx")
τy = FieldTimeSeries("surface_fluxes.jld2", "τy")
PE = FieldTimeSeries("surface_fluxes.jld2", "PE")

Qn  = @lift(interior(Q[$n],  :, :, 1))
τxn = @lift(interior(τx[$n], :, :, 1))
τyn = @lift(interior(τy[$n], :, :, 1))
PEn = @lift(interior(PE[$n], :, :, 1))

fig  = Figure(size = (1800, 800))
axQ  = Axis(fig[1, 1], title="Net heat flux Wm⁻²")
axPE = Axis(fig[1, 2], title="Net salt flux psu m s⁻¹")
axτx = Axis(fig[2, 1], title="Zonal wind stress Nm⁻²")
axτy = Axis(fig[2, 2], title="Meridional wind stress Nm⁻²")

hmQ  = heatmap!(axQ,  Qn,  colormap=:magma,   colorrange=(-800,  800))
hmPE = heatmap!(axPE, PEn, colormap=:haline,  colorrange=(-1e-5, 5e-5))
hmτx = heatmap!(axτx, τxn, colormap=:balance, colorrange=(-5e-4, 5e-4))
hmτy = heatmap!(axτy, τyn, colormap=:balance, colorrange=(-5e-4, 5e-4))

CairoMakie.record(fig, "surface_fluxes.mp4", 1:length(Q.times); framerate=5) do i 
    @info "doing $i of $(length(Q.times))"
    n[] = i
end

# Let's visualize the internal structure of temperature and salinity
x, y, z = nodes(ocean.model.tracers.T)

fig = Figure(size = (1200, 400))
axT = Axis(fig[1, 1], title="Internal temperature structure ᵒC")
axS = Axis(fig[1, 2], title="Internal salinity structure psu")

contourf!(axT, 1:48, z, interior(mean(ocean.model.tracers.T, dims=1), 1, :, :), colormap=:magma)
contourf!(axS, 1:48, z, interior(mean(ocean.model.tracers.S, dims=1), 1, :, :), colormap=:haline)

# # Running a high-resolution simulation
#
# What are the steps to modify the above script to run an eddying (quarter degree) simulation?
