module MRIOperatorsGPUArraysExt

using MRIOperators, GPUArrays, KernelAbstractions
using MRIOperators.FFTW, MRIOperators.LinearAlgebra
using MRIOperators.NFFT

include("ExplicitOp.jl")
include("Shutter.jl")
include("SensitivityOp.jl")
include("FieldmapNFFTOp.jl")

nfftParams(::NFFT.NFFTBackend, ::Type{<:AbstractGPUArray}) = (;)

end # module
