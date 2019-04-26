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

  # direct reco
  if recoParams[:reco] == "direct"
    shape, weights, cmap = setupDirectReco(acqData, recoParams)
    return reconstruction_direct_2d(acqData, shape, weights, cmap)
  end

  # iterative reco
  par = setupIterativeReco(acqData, recoParams)
  if recoParams[:reco] == "standard"
    return reconstruction_simple(acqData, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiEcho"
    return reconstruction_multiEcho(acqData, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoil"
    return reconstruction_multiCoil(acqData, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoilMultiEcho"
    return reconstruction_multiCoilMultiEcho(acqData, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  else
    @error "RecoModel $(recoParams[:reco]) not found."
  end

  return reconstruction_direct_2d(acqData, recoParams)
end

function reconstruction_3d(acqData::AcquisitionData, recoParams::Dict)
  recoParams = merge(defaultRecoParams(), recoParams)
  if recoParams[:reco] == "direct"
    shape, weights, cmap = setupDirectReco(acqData, recoParams)
    return reconstruction_direct_3d(acqData, shape, weights, cmap)
  end

  acqData2d = convert3dTo2d(acqData)
  par = setupIterativeReco(acqData2d, recoParams)
  if recoParams[:reco] == "standard"
    Ireco = reconstruction_simple(acqData2d, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiEcho"
    Ireco = reconstruction_multiEcho(acqData2d, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoil"
    Ireco = reconstruction_multiCoil(acqData2d, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoilMultiEcho"
    Ireco = reconstruction_multiCoilMultiEcho(acqData2d, par.shape, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  else
    @error "RecoModel $(recoParams[:reco]) not found."
  end
  Ireco = permutedims(Ireco,[3,1,2,4,5])

  return Ireco
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
