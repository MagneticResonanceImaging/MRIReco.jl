export reconstruction

################################################################################
include("DirectReconstruction.jl")
include("IterativeReconstruction.jl")
include("RecoParameters.jl")


"""
  reconstruction(acqData::AcquisitionData, recoParams::Dict)

The reconstruction method takes an AcquisitionData object and a parameter
dictionary and calculates an image from the given raw data.
"""
function reconstruction(acqData::AcquisitionData, recoParams::Dict)
  encodingDims = dims(trajectory(acqData))
  if encodingDims == 2 #"2D"
    return reconstruction_2d(acqData,recoParams)
  elseif encodingDims == 3 #"3D"
    return reconstruction_3d(acqData,recoParams)
  else
    error("$(encodingDims)d encoding not yet supported")
  end
  return reconstruction_2d(acqData,recoParams)
end


function reconstruction_2d(acqData::AcquisitionData, recoParams::Dict)
  recoParams = merge(defaultRecoParams(), recoParams)
  if recoParams[:reco] == "direct"
    return reconstruction_direct_2d(acqData, recoParams)
  elseif recoParams[:reco] == "standard"
    return reconstruction_simple(acqData, recoParams)
  elseif recoParams[:reco] == "multiEcho"
    return reconstruction_multiEcho(acqData, recoParams)
  elseif recoParams[:reco] == "multiCoil"
    return reconstruction_multiCoil(acqData, recoParams)
  elseif recoParams[:reco] == "multiCoilMultiEcho"
    return reconstruction_multiCoilMultiEcho(acqData, recoParams)
  else
    error("RecoModel $(recoParams[:reco]) not found.")
  end
  return reconstruction_direct_2d(acqData, recoParams)
end

function reconstruction_3d(acqData::AcquisitionData, recoParams::Dict)
  recoParams = merge(defaultRecoParams(), recoParams)
  if recoParams[:reco] == "direct"
    return reconstruction_direct_3d(acqData, recoParams)
  else
    error("3D reconstruction is not yet implimented")
  end
  return reconstruction_direct_3d(acqData, recoParams)
end

# This version stores the reconstructed data into a file
function reconstruction(acqData::AcquisitionData, recoParams::Dict, filename::String;
                        force=false)
  if !force && isfile(filename)
    return recoImage( RecoFileIBI(filename) )
  else
    I = reconstruction(acqData, recoParams)
    saveasRecoFile(filename, I, recoParams)
    return I
  end
end
