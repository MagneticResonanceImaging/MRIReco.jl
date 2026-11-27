# # [Undersample k-Space](@id undersample-kspace)
#
# Reconstruction methods for undersampled data are often tested on fully sampled data from which samples are removed
# afterwards. MRISampling does this for Cartesian data with
# [`sample_kspace`](@ref sample_kspace(::AcquisitionData, ::Any, ::AbstractString)). This guide shows how to choose a
# sampling pattern, keep the center of k-space, make the pattern reproducible and undersample 3D data and data with
# several contrasts. [Data Model](@ref data-model) explains how `AcquisitionData` stores the sampling pattern.
#
# ## Choosing a Pattern
#
# `sample_kspace` takes the data, the reduction factor and the name of a pattern and returns the undersampled data as
# a new object. The example keeps a quarter of the knee data from [Direct Reconstruction](@ref direct-reconstruction):

using MRIReco, MRIFiles, MRISampling

using Pkg.Artifacts
artifacts = joinpath(pkgdir(MRIReco), "docs", "Artifacts.toml")
datadir = ensure_artifact_installed("MRIRecoDocsData", artifacts; io = devnull)
acqData = AcquisitionData(ISMRMRDFile(joinpath(datadir, "knee_3dFSE_slice170.h5")))

acqDataSub = sample_kspace(acqData, 4.0, "poisson", seed = 1234)
size(acqData.kdata[1]), size(acqDataSub.kdata[1])

# By default, the pattern keeps complete profiles, here 80 of 320, as an accelerated acquisition would. The number of
# kept profiles is the number of profiles divided by the reduction factor, rounded down. Pass the factor as a
# floating-point number, such as `4.0`, because the random patterns do not accept an integer. There are six
# patterns:
#
# - `"regular"` keeps every n-th profile, as in a conventional parallel imaging acquisition. Use a whole number as
#   the factor, because other factors leave the end of k-space unsampled.
# - `"random"` chooses the profiles at random.
# - `"poisson"` chooses the profiles at random but keeps a minimum distance between them as long as there is room
#   left, which is known as Poisson-disk sampling.
# - `"vdPoissonGrid"` is a Poisson-disk pattern with variable density, which keeps more profiles near the center of
#   k-space. It is fast, also for large grids, and can keep a fully sampled center.
# - `"vdPoisson"` and `"vdPoissonBART"` are other Poisson-disk patterns with variable density: `"vdPoisson"` becomes
#   slow for large grids and cannot keep a fully sampled center, and `"vdPoissonBART"` uses the density law of the
#   BART toolbox. The variable-density patterns are described in more detail below.
#
# The docstring also lists `"lines"`, the pattern that `sample_kspace` uses internally to keep complete profiles. The
# following figure shows the first four patterns for the knee data, with the 30 central profiles always kept
# (`calsize`, see below; `"regular"` ignores it):

using CairoMakie

fig = Figure()
for (i, name) in enumerate(["regular", "random", "poisson", "vdPoissonGrid"])
    data = sample_kspace(acqData, 4.0, name, calsize = 30, seed = 1234)
    mask = zeros(320, 320)
    mask[data.subsampleIndices[1]] .= 1
    ax = Axis(fig[1, i], title = name, width = 200, height = 200,
              xlabel = "Sample in the readout", ylabel = "Phase encoding step")
    heatmap!(ax, mask, colormap = :grays)
end
resize_to_layout!(fig)
fig

# ## Keeping the Center of k-Space
#
# Some methods need a fully sampled region in the center of k-space. In
# [Multi-Coil Reconstruction](@ref multicoil-reconstruction), for example, ESPIRiT estimates the coil sensitivities
# from this region. For `"random"`, `"poisson"`, `"vdPoissonGrid"` and `"vdPoissonBART"`, the keyword `calsize` sets
# the number of central profiles that are always kept. The other profiles follow the pattern, and the reduction factor
# stays the same:

acqDataCal = sample_kspace(acqData, 4.0, "poisson", calsize = 30, seed = 1234)
profiles = unique(fld1.(acqDataCal.subsampleIndices[1], 320))
length(profiles), all(in(profiles), 146:175)

# `calsize` must not exceed the number of kept profiles. `"regular"` and `"vdPoisson"` ignore the keyword.
#
# ## Reproducible Patterns
#
# The random patterns differ from call to call. The keyword `seed` fixes them, so that a result can be reproduced:

a = sample_kspace(acqData, 4.0, "poisson", seed = 1234)
b = sample_kspace(acqData, 4.0, "poisson", seed = 1234)
a.subsampleIndices == b.subsampleIndices

# `seed` also resets Julia's default random number generator, so random numbers that are drawn afterwards repeat as
# well. To avoid this, pass a generator of your own as `rng`, which `seed` then seeds instead.
#
# ## 3D Data and Single Samples
#
# For 3D data, `sample_kspace` chooses complete readouts in the plane of the two phase encoding directions, so the
# patterns are two-dimensional. `"regular"` keeps every n-th point of this plane, counting along the first phase
# encoding direction. `sample` returns such a pattern without applying it to data, as a vector with the indices of
# the kept points of a grid:

idx = sample((64, 64), 4.0, "poisson", calsize = 12, seed = 1234)
length(idx)

# In two dimensions, `calsize` keeps a central square, here of 12 × 12 points, which must not contain more points
# than the pattern keeps. The following figure shows a Poisson-disk pattern with uniform and with variable density:

fig = Figure()
for (i, (name, options)) in enumerate(["poisson" => (calsize = 12,), "vdPoissonGrid" => (calsize = 12,)])
    mask = zeros(64, 64)
    mask[sample((64, 64), 4.0, name; seed = 1234, options...)] .= 1
    ax = Axis(fig[1, i], title = name, width = 250, height = 250,
              xlabel = "Phase encoding step 1", ylabel = "Phase encoding step 2")
    heatmap!(ax, mask, colormap = :grays)
end
resize_to_layout!(fig)
fig

# With `profiles = false`, `sample_kspace` uses such a two-dimensional pattern for 2D data as well and removes single
# samples instead of complete profiles, which treats a slice like the phase encoding plane of a 3D acquisition. The
# random patterns support this option only for 2D data.
#
# ## Variable-Density Poisson-Disk Patterns
#
# `"vdPoissonGrid"` and `"vdPoissonBART"` are Poisson-disk patterns whose density decreases from the center to the
# edge of k-space. Unlike `"vdPoisson"`, they keep a fully sampled center with `calsize` and remain fast for large
# grids. Both keep exactly the number of samples given by the reduction factor:
#
# - `"vdPoissonGrid"` is the pattern to start with. It places the samples directly on the grid, with a minimum
#   distance that grows linearly with the distance to the center at a rate set by the keyword `vardens` (2 by
#   default; 0 gives a uniform pattern).
# - `"vdPoissonBART"` uses the density law of the `poisson` command of [BART](https://github.com/mrirecon/bart), for
#   example to compare with patterns made with BART: the minimum distance between the samples grows with the square
#   of the distance to the center, and the keyword `vardens` sets how fast (20 by default, as the option `-v` of
#   BART; 0 gives a uniform pattern).
#
# A strong variable density does not allow small reduction factors, because the edge of k-space cannot be denser than
# the grid. Both patterns then throw an error with the smallest reachable factor, for example about 8 for a
# two-dimensional `"vdPoissonBART"` pattern with `vardens = 20` and `elliptical = true`; use a smaller `vardens` in
# this case. Both patterns call the strength of the variable density `vardens`, but on different scales, because
# their density laws differ.
#
# ### Profiles of 2D Data
#
# With 2D data and complete profiles, the patterns are one-dimensional: they choose the kept phase encoding steps, with
# more steps near the center of k-space and a fully sampled central block of `calsize` steps. The example keeps a
# quarter of the profiles of the knee data, including the 30 central ones:

fig = Figure()
for (i, (name, options)) in enumerate(["vdPoissonGrid" => (vardens = 2.0,), "vdPoissonBART" => (vardens = 20.0,)])
    data = sample_kspace(acqData, 4.0, name; calsize = 30, seed = 1234, options...)
    mask = zeros(320, 320)
    mask[data.subsampleIndices[1]] .= 1
    ax = Axis(fig[1, i], title = name, width = 250, height = 250,
              xlabel = "Sample in the readout", ylabel = "Phase encoding step")
    heatmap!(ax, mask, colormap = :grays)
end
resize_to_layout!(fig)
fig

# ### Phase Encoding Plane of 3D Data
#
# For 3D data, the patterns are two-dimensional and sample the plane of the two phase encoding directions. `sample`
# returns such a pattern for a 256 × 256 plane, here with a reduction factor of 8, a central square of 24 × 24
# samples, and `elliptical = true`, which leaves the corners of the plane empty:

fig = Figure()
for (i, (name, options)) in enumerate(["vdPoissonGrid" => (vardens = 2.0,), "vdPoissonBART" => (vardens = 20.0,)])
    mask = zeros(256, 256)
    mask[sample((256, 256), 8.0, name; calsize = 24, elliptical = true, seed = 1234, options...)] .= 1
    ax = Axis(fig[1, i], title = name, width = 250, height = 250,
              xlabel = "Phase encoding step 1", ylabel = "Phase encoding step 2")
    heatmap!(ax, mask, colormap = :grays)
end
resize_to_layout!(fig)
fig

# [`sample_lines`](@ref sample_lines(::Tuple{Int64, Int64, Int64}, ::Any)) applies such a pattern to all samples of
# the readout and returns the indices in the complete 3D grid, here of 16 readout samples and a 64 × 64 plane:

idx = sample_lines((16, 64, 64), 8.0; sampleFunc = "vdPoissonGrid", calsize = 12, elliptical = true, seed = 1234)
length(idx) == 16 * floor(Int, 64^2 / 8.0)

#
# ## Several Contrasts
#
# For data with several contrasts, each contrast gets its own random pattern by default, and `rand = false` uses one
# pattern for all contrasts. [Multi-Echo Reconstruction](@ref multiecho-reconstruction) shows how this choice affects
# a joint reconstruction of the contrasts.
#
# ## Other Kinds of Data
#
# `sample_kspace` expects fully sampled Cartesian data whose samples are stored in the order of the k-space grid. For
# fully sampled data, this is the case if the indices of all contrasts are sorted, as for the knee data:

all(issorted, acqData.subsampleIndices)

# - Data that is loaded from a file can be stored in another order. `sample_kspace` then assigns the kept samples to
#   wrong k-space positions without an error
#   ([#352](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/352)). Bring such data into the order of
#   the grid first with `AcquisitionData(kDataCart(acqData))`, and add `enc2D = true` for data with several slices.
#   This resets the field of view and the sampling times, and it treats samples that are exactly zero in the first
#   coil of the first slice and repetition as not measured.
# - Data that is already undersampled cannot be undersampled again with `sample_kspace`.
# - For non-Cartesian data, `sample_kspace` removes complete profiles from the k-space data but keeps all nodes of
#   the trajectory ([#353](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/353)). Pass the result to
#   [`convertUndersampledData`](@ref), which reduces the trajectory to the nodes of the kept samples:
#   `convertUndersampledData(sample_kspace(acqData, 4.0, "random"))`.
# - To apply a pattern of your own, set the samples that are not measured to zero in the k-space array, as shown in
#   [Extract Cartesian k-Space](@ref extract-kspace).
