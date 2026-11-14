module MRIRecoAMDGPUExt

using MRIReco, AMDGPU

MRIReco.set_device!(::Type{<:ROCArray}, index) = AMDGPU.device!(index) 

end