export reconstruction

################################################################################
include("DirectReconstruction.jl")
include("IterativeReconstruction.jl")
include("RecoParameters.jl")


"""
    reconstruction(acqData::AcquisitionData, recoParams::Dict)

takes an AcquisitionData object and a parameter
dictionary and calculates an image from the given raw data.

Reconstruction types are specified by the symbol `:reco`.
Valid reconstruction names are:
* :direct - direct Fourier reconstruction
* :standard           - iterative reconstruction for all contrasts, coils & slices independently
* :multiEcho          - iterative joint reconstruction of all echo images
* :multiCoil          - SENSE-type iterative reconstruction
* :multiCoilMultiEcho - SENSE-type iterative reconstruction of all echo images
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

"""
    reconstruction_2d(acqData::AcquisitionData, recoParams::Dict)

Performs image reconstruction of a 2d encoded AcquisitionData object.
Parameters are specified in a dictionary.

Reconstruction types are specified by the symbol `:reco`.
Valid reconstruction names are:
* :direct - direct Fourier reconstruction
* :standard           - iterative reconstruction for all contrasts, coils & slices independently
* :multiEcho          - iterative joint reconstruction of all echo images
* :multiCoil          - SENSE-type iterative reconstruction
* :multiCoilMultiEcho - SENSE-type iterative reconstruction of all echo images
"""
function reconstruction_2d(acqData::AcquisitionData, recoParams::Dict)
  recoParams = merge(defaultRecoParams(), recoParams)

  # direct reco
  if recoParams[:reco] == "direct"
    reconSize, weights, cmap = setupDirectReco(acqData, recoParams)
    return reconstruction_direct_2d(acqData, reconSize, weights, cmap)
  end

  # iterative reco
  par = setupIterativeReco(acqData, recoParams)
  if recoParams[:reco] == "standard"
    return reconstruction_simple(acqData, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiEcho"
    return reconstruction_multiEcho(acqData, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoil"
    return reconstruction_multiCoil(acqData, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoilMultiEcho"
    return reconstruction_multiCoilMultiEcho(acqData, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  else
    @error "RecoModel $(recoParams[:reco]) not found."
  end

  return reconstruction_direct_2d(acqData, recoParams)
end

"""
    reconstruction_3d(acqData::AcquisitionData, recoParams::Dict)

Performs image reconstruction of a 3d encoded AcquisitionData object.
Parameters are specified in a dictionary.

Reconstruction types are specified by the symbol `:reco`.
Valid reconstruction names are:
* :direct - direct Fourier reconstruction
* :standard           - iterative reconstruction for all contrasts, coils & slices independently
* :multiEcho          - iterative joint reconstruction of all echo images
* :multiCoil          - SENSE-type iterative reconstruction
* :multiCoilMultiEcho - SENSE-type iterative reconstruction of all echo images
"""
function reconstruction_3d(acqData::AcquisitionData, recoParams::Dict)
  recoParams = merge(defaultRecoParams(), recoParams)
  if recoParams[:reco] == "direct"
    reconSize, weights, cmap = setupDirectReco(acqData, recoParams)
    return reconstruction_direct_3d(acqData, reconSize, weights, cmap)
  end

  acqData2d = convert3dTo2d(acqData)
  par = setupIterativeReco(acqData2d, recoParams)
  if recoParams[:reco] == "standard"
    Ireco = reconstruction_simple(acqData2d, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiEcho"
    Ireco = reconstruction_multiEcho(acqData2d, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoil"
    Ireco = reconstruction_multiCoil(acqData2d, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  elseif recoParams[:reco] == "multiCoilMultiEcho"
    Ireco = reconstruction_multiCoilMultiEcho(acqData2d, par.reconSize, par.reg, par.sparseTrafo, par.weights, par.solvername, par.senseMaps, par.correctionMap, par.method, par.normalize, recoParams)
  else
    @error "RecoModel $(recoParams[:reco]) not found."
  end
  Ireco = permutedims(Ireco,[3,1,2,4,5])

  return Ireco
end

"""
    reconstruction(acqData::AcquisitionData, recoParams::Dict,filename::String; force=false)

performs the same image reconstrucion as `reconstruction(acqData::AcquisitionData, recoParams::Dict)`
and saves the image in a file with name `filename`.
If `force=false`, the reconstructed image is loaded from the the file `filename` if the latter is
present.
"""
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
