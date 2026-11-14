module MRIRecoCUDAExt

using MRIReco, CUDA

MRIReco.set_device!(::Type{<:CuArray}, index) = CUDA.device!(index) 

end