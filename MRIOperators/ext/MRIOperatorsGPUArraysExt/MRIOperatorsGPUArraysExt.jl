module MRIOperatorsGPUArraysExt

using MRIOperators, GPUArrays, MRIOperators.FFTW

include("ExplicitOp.jl")
include("Shutter.jl")
include("SensitivityOp.jl")


end # module
