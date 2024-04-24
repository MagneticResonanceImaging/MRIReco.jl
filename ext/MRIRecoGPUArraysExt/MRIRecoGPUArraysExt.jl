module MRIRecoGPUArraysExt

using MRIReco, GPUArrays, MRIReco.FLoops

MRIReco.executor(::Type{<:AbstractGPUArray}) = SequentialEx()
MRIReco.copyOpsFn(::Type{<:AbstractGPUArray}) = identity

end