# Getting Started

## Installation

MRIReco.jl requires Julia 1.10 or newer. Install the main package with

```julia
using Pkg
Pkg.add("MRIReco")
```

`using MRIReco` also loads MRIBase, MRIOperators and [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/), which provides the solvers. The other packages are installed in the same way when you need them, for example `Pkg.add("MRIFiles")` to read raw data files. The examples in this documentation also use [CairoMakie](https://docs.makie.org/stable/) for plots and [ImageUtils](https://github.com/IBIResearch/ImageUtils.jl) for test images.

### Development Version

The [development documentation](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/) describes the `master` branch, which can contain changes that are not released yet. To use `master`, add the package and the packages of this repository that it depends on, each with `rev = "master"`. Otherwise Pkg takes these dependencies from the registry, and their released versions may not work with `master`. For MRIReco, which depends on MRIBase and MRIOperators:

```julia
using Pkg
Pkg.add([PackageSpec(name = "MRIBase", rev = "master"),
         PackageSpec(name = "MRIOperators", rev = "master"),
         PackageSpec(name = "MRIReco", rev = "master")])
```

For other packages, such as MRIFiles or MRISimulation, the [table of packages](@ref "Packages") on the Home page lists the packages of this repository that they depend on. `Pkg.update()` moves the packages to the newest commit of `master`, and `Pkg.free(["MRIBase", "MRIOperators", "MRIReco"])` returns them to their released versions. To work on the packages themselves, follow the [development setup](@ref "Development Setup") of the contributing guide.

## A First Reconstruction

This example simulates k-space data of a Shepp-Logan phantom along a radial trajectory and reconstructs an image from it. Besides MRIReco it needs MRISimulation, ImageUtils and CairoMakie.

```@example gettingstarted
using MRIReco, MRISimulation, CairoMakie
using ImageUtils: shepp_logan

N = 256
phantom = shepp_logan(N)

params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Radial"
params[:numProfiles] = floor(Int64, pi/2 * N)
params[:numSamplingPerProfile] = 2N
acqData = simulation(phantom, params)
nothing # hide
```

`acqData` is an [`AcquisitionData`](@ref) object that holds the simulated k-space data together with the trajectory. A direct reconstruction applies the adjoint of the Fourier encoding to the density-compensated data:

```@example gettingstarted
params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N, N)
image = reconstruction(acqData, params)
size(image)
```

The reconstructed image always has six dimensions: x, y, z, echoes, coils and repetitions. Here only the first two have more than one entry:

```@example gettingstarted
images = [phantom, abs.(image[:, :, 1, 1, 1, 1])]
titles = ["Phantom", "Reconstruction"]
fig = Figure(size = (700, 350))
for i in 1:2
    ax = Axis(fig[1, i], aspect = DataAspect(), yreversed = true, title = titles[i])
    heatmap!(ax, images[i]', colormap = :grays)
    hidedecorations!(ax)
end
fig
```

## Next Steps

- The tutorials walk through complete reconstructions step by step, starting with [Direct Reconstruction](@ref direct-reconstruction) and [Iterative Reconstruction](@ref iterative-reconstruction).
- The how-to guides show how to solve specific tasks, such as binning profiles, extracting Cartesian k-space data or [running a reconstruction on the GPU](@ref gpu-reconstruction).
- The API reference lists the functions and types of each package.
