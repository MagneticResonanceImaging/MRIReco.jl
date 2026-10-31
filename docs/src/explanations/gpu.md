# [GPU Acceleration](@id gpu-acceleration)

An iterative reconstruction spends most of its time in applying the encoding operator and its adjoint, which consist of Fourier transforms and multiplications with coil sensitivities, and in the vector operations of the solver. These operations act on all voxels or samples at once and therefore suit the many cores of a GPU. MRIReco.jl runs reconstructions on GPUs without GPU-specific code in the reconstruction functions. This page explains how, and [Run a Reconstruction on the GPU](@ref gpu-reconstruction) shows how to use it. The approach is presented in more detail in the talk [Generic GPU-Acceleration for Medical Image Reconstruction](https://www.youtube.com/watch?v=SsIwveAJOJ4) at JuliaCon 2026, whose [slides](https://pretalx.com/media/juliacon-2026/submissions/F8GRBR/resources/JuliaCon2026_dYxDkif.pdf) are available as well.

## Generic Array Types

Julia's GPU packages, for example [CUDA.jl](https://cuda.juliagpu.org/stable/) for NVIDIA GPUs and [AMDGPU.jl](https://amdgpu.juliagpu.org/stable/) for AMD GPUs, provide array types whose data is stored on the GPU, such as `CuArray` and `ROCArray`. Broadcasting, reductions, linear algebra and FFTs run on the GPU when they are called with these arrays. Code that works with any `AbstractArray` and doesn't assume that its arrays are `Array`s therefore runs on the GPU as well.

The packages that MRIReco.jl [builds on](architecture.md) are written this way: the solvers and regularization terms of [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/), the general operators of [LinearOperatorCollection.jl](https://juliaimagerecon.github.io/LinearOperatorCollection.jl/stable/) and the encoding operators of MRIOperators. None of them depends on a GPU package.

## The State Follows the Data

The reconstruction parameter `arrayType` decides where a reconstruction runs. MRIReco converts the measured data, the density compensation weights and the coil sensitivities to this array type. From the array type and the precision of the data, it also derives the storage type of the operators, for example a `CuArray` with element type `ComplexF32` for single-precision data. This is why `arrayType` is given without parameters. The operators allocate their internal arrays with this storage type, so these arrays are on the GPU as well. The solvers create their state, such as the current estimate and the residual, with `similar` from the measured data, so it is stored where the data is.

The rest runs on the CPU: reading the data and estimating the coil sensitivities before the reconstruction, and computing the density compensation weights at its start, before MRIReco converts them. At the end, MRIReco copies the image back to the CPU and stores it in an `AxisArray`, as for a reconstruction on the CPU.

## GPU Code in Package Extensions

Some operations can't be written efficiently with generic array operations. For these, MRIOperators contains GPU kernels, which are written with [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) and are therefore not specific to one GPU vendor: the multiplication with coil sensitivities in `SensitivityOp`, the off-resonance correction in `FieldmapNFFTOp`, the explicit encoding operator `ExplicitOp` and the circular shutter that is applied for trajectories with circular k-space coverage, such as radial and spiral trajectories. Non-Cartesian Fourier transforms use the GPU implementation of [NFFT.jl](https://juliamath.github.io/NFFT.jl/stable/).

This code is in [package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)): the kernels in `MRIOperatorsGPUArraysExt` of MRIOperators, the GPU settings of MRIReco that are described below in `MRIRecoGPUArraysExt`, and the GPU implementation of the NFFT in an extension of NFFT.jl. These extensions depend on [GPUArrays.jl](https://github.com/JuliaGPU/GPUArrays.jl), KernelAbstractions.jl and Adapt.jl, which the GPU packages build on, and Julia loads them as soon as these packages are loaded. Loading CUDA.jl therefore loads the extensions as well, and without a GPU package they are not loaded at all.

## Parallelism on the CPU and the GPU

The [MRIReco.jl paper](https://doi.org/10.1002/mrm.28792) describes multi-threading on two levels. Within an operator, for example the multi-coil encoding operator, several Fourier transforms are computed in parallel. Above that, the reconstruction types `"multiEcho"`, `"multiCoil"`, `"multiCoilMultiEcho"` and `"multiCoilMultiEchoSubspace"` divide the data into independent problems, for example one per slice and repetition, which can be solved in parallel. On the CPU, MRIReco copies the Fourier operator for each coil, so that the copies can be applied in parallel, and solves the independent problems with the [OhMyThreads.jl](https://juliafolds2.github.io/OhMyThreads.jl/stable/) scheduler `DynamicScheduler`.

On the GPU, each operation of a single problem already runs in parallel. MRIReco therefore solves the independent problems one after another, with the `SerialScheduler`, and doesn't copy the Fourier operator, which saves GPU memory. The parameter `scheduler` changes the scheduler for the independent problems. With several Julia threads, the shared Fourier operator currently gives wrong results ([#317](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/317)).

## Further Reading

- The GPU guides of [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/generated/howto/gpu_acceleration/) and [LinearOperatorCollection.jl](https://juliaimagerecon.github.io/LinearOperatorCollection.jl/stable/generated/howtos/gpu/) show how the solvers and operators are used on the GPU directly.
- [JuliaGPU](https://juliagpu.org) gives an overview of the GPU packages of Julia.
