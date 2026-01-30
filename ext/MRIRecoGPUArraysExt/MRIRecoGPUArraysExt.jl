module MRIRecoGPUArraysExt

using MRIReco, GPUArrays, MRIReco.OhMyThreads

MRIReco.scheduler(::Type{<:AbstractGPUArray}) = SerialScheduler()
MRIReco.copyOpsFn(::Type{<:AbstractGPUArray}) = identity

end