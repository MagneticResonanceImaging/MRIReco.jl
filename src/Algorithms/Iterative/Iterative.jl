setupIterativeReconSize(acqData::AcquisitionData, ::Nothing) = setupIterativeReconSize(acqData, encodingSize(acqData))
function setupIterativeReconSize(acqData::AcquisitionData, reconSize::NTuple)
  encodingDims = ndims(trajectory(acqData))
  if encodingDims==3 && numSlices(acqData)>1
    @error "reconstruction of multiple 3d-encoded volumina is not yet supported"
  end

  red3d = ndims(trajectory(acqData,1))==2 && length(reconSize)==3
  if red3d  # acqData is 3d data converted to 2d
    reconSize = (reconSize[2], reconSize[3])
  end
  return reconSize
end


include("Standard.jl")
include("MultiEcho.jl")
include("MultiCoil.jl")
include("MultiCoilMultiEcho.jl")
include("MultiCoilMultiEchoSubspace.jl")
