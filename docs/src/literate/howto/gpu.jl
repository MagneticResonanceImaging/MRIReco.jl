# # [Run a Reconstruction on the GPU](@id gpu-reconstruction)
#
# MRIReco.jl can run reconstructions on a GPU. This guide shows how to do this on an NVIDIA GPU with
# [CUDA.jl](https://cuda.juliagpu.org/stable/). Whether the GPU is faster than the CPU depends on the size of the
# problem and on the hardware. Iterative reconstructions, which apply the encoding operator many times, benefit the
# most. [GPU Acceleration](@ref gpu-acceleration) explains how MRIReco.jl runs on the GPU.
#
# ## Installing CUDA.jl
#
# MRIReco.jl loads its GPU code together with a GPU package, which has to be installed in the same environment as
# MRIReco.jl. For NVIDIA GPUs, this is CUDA.jl:
# ```julia
# using Pkg
# Pkg.add("CUDA")
# ```
# CUDA.jl only needs the driver of the GPU and downloads the CUDA libraries itself. After `using CUDA`,
# `CUDA.functional()` returns `true` if CUDA.jl can use the GPU, and `CUDA.versioninfo()` lists the GPUs and the
# versions of the libraries.
#
# ## Choosing the Array Type
#
# The parameter `arrayType` sets the type of the arrays that the reconstruction computes with. The default `Array`
# keeps the reconstruction on the CPU, and the GPU array type `CuArray` moves it to the GPU:
# ```julia
# using MRIReco, CUDA
# params[:arrayType] = CuArray
# ```
# Pass the type without parameters. MRIReco.jl adds the element type, which follows from the data, so for example
# `CuArray{ComplexF32}` fails.
#
# Apart from this parameter, the reconstruction is set up as on the CPU. The documentation is built on a computer
# without a GPU, so the following example uses `gpu = Array` and runs on the CPU. With an NVIDIA GPU,
# `gpu = CuArray` runs it on the GPU. The example reconstructs the knee data from
# [Direct Reconstruction](@ref direct-reconstruction) with CG-SENSE:

using MRIReco, MRIFiles, MRICoilSensitivities
gpu = Array # with an NVIDIA GPU: using CUDA; gpu = CuArray

datadir = joinpath(pkgdir(MRIReco), "docs", "src", "examples", "data")
acqData = AcquisitionData(ISMRMRDFile(joinpath(datadir, "knee_3dFSE_slice170.h5")))

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (320, 320)
params[:senseMaps] = espirit(acqData)
params[:solver] = CGNR
params[:reg] = L2Regularization(1e-3)
params[:normalizeReg] = MeasurementBasedNormalization()
params[:iterations] = 10
params[:arrayType] = gpu
image = reconstruction(acqData, params)
nothing #hide

# Only the reconstruction itself runs on the GPU. Reading the data and estimating the coil sensitivities with
# [`espirit`](@ref) happen on the CPU, and `reconstruction` copies the image back to the CPU. It returns the same
# six-dimensional `AxisArray` as without GPU, so the code after the reconstruction stays the same. Because the GPU
# computes in a different order, its result differs from the CPU result by rounding errors. With non-Cartesian data,
# the difference is larger, because the NFFT evaluates its window function differently on the CPU and the GPU, but it
# is still far smaller than the error of the reconstruction itself.
#
# ## Precision
#
# The reconstruction computes in the precision of the data. MRIFiles reads raw data in single precision, as the type
# parameter of the knee data shows:

typeof(acqData)

# Most GPUs compute much faster in single than in double precision, so single-precision data is preferable on the
# GPU.
#
# ## Selecting a GPU
#
# CUDA.jl uses the first GPU by default. On a computer with several GPUs, select another one with `CUDA.device!`
# before the reconstruction, for example the second GPU with `CUDA.device!(1)`.
#
# ## Current Limitations
#
# Some parts of MRIReco.jl don't work on the GPU yet:
#
# - With several Julia threads, `"multiCoil"` and `"multiCoilMultiEcho"` give wrong results on the GPU
#   ([#317](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/317)). Start Julia with one thread, which is
#   the default, for these reconstructions.
# - Noise decorrelation with `:noiseData` fails on the GPU with Julia 1.12 and 1.13, and for `"multiCoilMultiEcho"`
#   also with Julia 1.10 ([#319](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/319)).
# - `"multiCoilMultiEchoSubspace"` fails on the GPU
#   ([#320](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/320)).
# - The GPU code is not specific to CUDA.jl, but with AMD GPUs and [AMDGPU.jl](https://amdgpu.juliagpu.org/stable/),
#   the iterative reconstructions currently fail
#   ([#318](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/318)).
