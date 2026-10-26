# # [Data Model](@id data-model)
#
# MRIReco.jl represents measured MRI data with a few data types, which are defined in MRIBase and shared by all
# packages of the framework. This page explains what these types contain, how data flows between them and how
# k-space data and images are laid out.
#
# ![Data flow: MRIFiles reads a raw data file into RawAcquisitionData and writes it back. RawAcquisitionData is
# converted into AcquisitionData and back. AcquisitionData can also be created by simulation with MRISimulation.
# The function reconstruction from MRIReco turns AcquisitionData and a dictionary of parameters into an image stored
# as an AxisArray.](../../assets/data-flow.svg)
#
# Raw data is read from a file into a [`RawAcquisitionData`](@ref) object and converted into an
# [`AcquisitionData`](@ref) object, which the [`reconstruction`](@ref) function turns into an image. `AcquisitionData`
# can also be created by [`simulation`](@ref). In the diagram, functions are drawn with round ends, and arrows with
# two heads mark conversions that work in both directions.
#
# The examples on this page use one slice of a knee dataset from [mridata.org](http://mridata.org), which is
# included in the MRIReco.jl repository.

using MRIReco, MRIFiles

datadir = joinpath(pkgdir(MRIReco), "docs", "src", "examples", "data")
filename = joinpath(datadir, "knee_3dFSE_slice170.h5")
nothing #hide

# ## Raw Data
#
# `RawAcquisitionData` closely follows the [ISMRMRD](https://ismrmrd.readthedocs.io/) format, an open standard for
# MRI raw data. It has two fields: `params`, a dictionary with the information from the XML header of an ISMRMRD
# file, and `profiles`, a vector with one [`Profile`](@ref) for each readout. MRIFiles reads ISMRMRD files and data
# from Bruker scanners into this type and writes it to ISMRMRD files. Data from other vendors can be converted into
# ISMRMRD files with the converters of the ISMRMRD project, such as
# [`siemens_to_ismrmrd`](https://github.com/ismrmrd/siemens_to_ismrmrd).

rawData = RawAcquisitionData(ISMRMRDFile(filename))
nothing #hide

# The XML header is nested, but `params` stores it flattened, so that each parameter can be accessed directly. The
# matrix size of the encoded k-space, for example, is stored in the header under `encoding/encodedSpace/matrixSize`
# and is available as

rawData.params["encodedSize"]

# Each profile holds the samples of one readout and, for non-Cartesian acquisitions, their k-space positions in its
# field `traj`. The samples of all coils are stored in its field `data`, with one row per sample and one column per
# coil:

profile = rawData.profiles[1]
size(profile.data)

# The header of a profile, an [`AcquisitionHeader`](@ref), describes the readout. Its encoding counters `idx`
# give the position of the readout in the acquisition, such as the phase encoding step, the slice, the contrast and
# the repetition. Like in ISMRMRD, the counters start at zero. The 320 profiles of the knee data cover the phase
# encoding steps 0 to 319:

extrema(Int(p.head.idx.kspace_encode_step_1) for p in rawData.profiles)

# The header also contains flags that mark special readouts, for example noise measurements or calibration data.
# The counters and flags are described in the
# [ISMRMRD documentation](https://ismrmrd.readthedocs.io/en/latest/mrd_raw_data.html).

# ## Preprocessed Data
#
# The profiles in a raw data file are not necessarily stored in the order of the k-space positions, and their layout
# depends on the sequence. For reconstruction, the data is therefore converted into an `AcquisitionData` object:

acqData = AcquisitionData(rawData)
nothing #hide

# The conversion groups the profiles by contrast, slice and repetition, averages repeated measurements, removes
# samples that the header marks to be discarded and reverses readouts that were acquired in reverse direction. It
# also creates the sampling trajectory and determines the trajectory node of each sample. The conversion works in
# the other direction as well: [`RawAcquisitionData(::AcquisitionData)`](@ref) creates raw data, for example to
# store simulated data in an ISMRMRD file.
#
# `AcquisitionData` is defined as
#
# ```julia
# mutable struct AcquisitionData{T <: AbstractFloat, D}
#   sequenceInfo::Dict{Symbol,Any}
#   traj::Vector{Trajectory{T}}
#   kdata::Array{Matrix{Complex{T}},3}
#   subsampleIndices::Vector{Vector{Int64}}
#   encodingSize::NTuple{D, Int64}
#   fov::NTuple{3,Float64}
# end
# ```
#
# where `T` is the floating point precision and `D` the number of encoding dimensions, 2 or 3:

typeof(acqData)

# The k-space data `kdata` is an array of matrices. Its three dimensions are the contrasts (or echoes), the slices
# and the repetitions. Each matrix holds the data of one contrast, slice and repetition, with one row per k-space
# sample and one column per coil:

size(acqData.kdata), size(acqData.kdata[1])

# The knee data has one contrast, slice and repetition, and 320 × 320 = 102400 samples from 8 coils. The functions
# [`numContrasts`](@ref MRIBase.numContrasts), [`numSlices`](@ref), [`numRepetitions`](@ref) and
# [`numChannels`](@ref numChannels(::AcquisitionData)) return these sizes.
#
# `sequenceInfo` can store additional information on the pulse sequence and is empty for the knee data.
# `encodingSize` is the size of the image grid that the data encodes, and the reconstruction uses it as the default
# image size. It is needed in particular for non-Cartesian trajectories, where the grid size cannot be derived from
# the trajectory. `fov` is the size of the encoded field of view in millimeters. It always has three entries, and for
# this knee slice, which was taken from a 3D acquisition, the third entry is the extent of the 3D volume in the third
# dimension:

acqData.encodingSize, acqData.fov

# ## Trajectories and Undersampling
#
# The field `traj` of `AcquisitionData` contains one [`Trajectory`](@ref) for each contrast. Its most important field
# is `nodes`, a matrix with the k-space positions of the samples, with one row per dimension and one column per node.
# The positions are measured in units of the inverse pixel size of the encoding grid, so that a fully sampled grid
# covers the interval `[-0.5, 0.5)` in each dimension. [`kspaceNodes`](@ref) returns this matrix:

tr = trajectory(acqData, 1)
size(kspaceNodes(tr)), extrema(kspaceNodes(tr))

# The knee data was acquired on a Cartesian grid. MRIBase also provides constructors for non-Cartesian trajectories,
# such as radial and spiral trajectories and, in 3D, kooshball and stack-of-stars trajectories. The following figure
# compares a Cartesian, a radial and a spiral trajectory:

using CairoMakie

fig = Figure(size = (900, 320))
for (i, (name, numProfiles, numSamples)) in
        enumerate([("Cartesian", 16, 16), ("Radial", 16, 32), ("Spiral", 4, 256)])
    nodes = kspaceNodes(trajectory(Float64, name, numProfiles, numSamples))
    ax = Axis(fig[1, i], title = name, aspect = DataAspect(),
              limits = (-0.5, 0.5, -0.5, 0.5), xlabel = "kx", ylabel = "ky")
    scatter!(ax, nodes[1, :], nodes[2, :], markersize = 4)
end
fig

# `subsampleIndices` contains, for each contrast, the indices of the trajectory nodes at which data was measured,
# in the order of the rows of `kdata`. For fully sampled data, these are all nodes. Undersampling, for example with
# [`sample_kspace`](@ref sample_kspace(::AcquisitionData, ::Any, ::AbstractString)) from MRISampling, keeps the
# trajectory and removes the samples that are not measured from `kdata` and `subsampleIndices`:

using MRISampling

acqDataUnder = sample_kspace(acqData, 2.0, "regular")
length(acqDataUnder.subsampleIndices[1]), size(acqDataUnder.kdata[1])

# Since `traj` and `subsampleIndices` have one entry per contrast, the trajectory and the undersampling pattern can
# differ between contrasts, but they are shared by all slices and repetitions.

# ## Images
#
# The reconstruction returns the image as an `AxisArray` from
# [AxisArrays.jl](https://github.com/JuliaArrays/AxisArrays.jl), an array whose dimensions have names and
# coordinates. It always has six dimensions. Their names and coordinates can be read with the functions of
# AxisArrays.jl, which has to be added to the environment to be loaded:

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (320, 320)
image = reconstruction(acqData, params)

using AxisArrays
axisnames(image)

# The first three dimensions are spatial, and for 2D data the slices are stored along `z`. The fourth dimension
# holds the contrasts, the fifth the coils and the sixth the repetitions. The direct reconstruction used here
# reconstructs each coil separately, so the coil dimension has length 8. Reconstructions that combine the coils with
# coil sensitivities return a single image in this dimension. The spatial axes have coordinates in millimeters,
# computed from `fov` and `encodingSize`:

axisvalues(image)[1]

# `parent(image)` returns the underlying array. MRIReco does not write image files itself; images can be stored
# with other packages, for example in the NIfTI format with `saveImage` from
# [ImageUtils.jl](https://github.com/IBIResearch/ImageUtils.jl), which keeps the pixel spacing.

# ## Array Layouts
#
# The order of the dimensions differs between the data types, which matters when data is passed between them:
#
# | Data | Dimensions |
# |:---|:---|
# | `kdata` of `AcquisitionData` | contrast, slice, repetition; each matrix: sample, coil |
# | k-space array of [`kDataCart`](@ref) | x, y, z or slice, coil, contrast, repetition |
# | reconstructed image | x, y, z or slice, contrast, coil, repetition |
# | coil sensitivities of [`espirit`](@ref) | x, y, z or slice, coil |
#
# A k-space array in the layout of `kDataCart` can also be converted into `AcquisitionData`. Note that the coil and
# contrast dimensions are in a different order in k-space arrays and in images:

size(kDataCart(acqData)), size(image)
