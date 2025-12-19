module MRIOperatorsGPUArraysExt

using MRIOperators, GPUArrays, KernelAbstractions
using MRIOperators.FFTW, MRIOperators.LinearAlgebra
using MRIOperators.NFFT
using MRIOperators.OhMyThreads

include("ExplicitOp.jl")
include("Shutter.jl")
include("SensitivityOp.jl")
include("FieldmapNFFTOp.jl")

# Flags for GPUs seem to break dispatch
MRIOperators.fftParams(::Type{<:AbstractGPUArray}) = (;)
MRIOperators.nfftParams(::NFFT.NFFTBackend, ::Type{<:AbstractGPUArray}) = (;)

end # module
