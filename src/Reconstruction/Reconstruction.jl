export reconstruction

include("RecoParameters.jl")
include("DirectReconstruction.jl")
include("IterativeReconstruction.jl")


"""
    reconstruction(acqData::AcquisitionData, recoParams::Dict)

Performs image reconstruction of an AcquisitionData object.
Parameters are specified in a dictionary.

Reconstruction types are specified by the symbol `:reco`.
Valid reconstruction names are:
* :direct - direct Fourier reconstruction
* :standard           - iterative reconstruction for all contrasts, coils & slices independently
* :multiEcho          - iterative joint reconstruction of all echo images
* :multiCoil          - SENSE-type iterative reconstruction
* :multiCoilMultiEcho - SENSE-type iterative reconstruction of all echo images
"""
function reconstruction(acqData::AcquisitionData, recoParams::Dict)
  # check dimensionality of encoding
  encodingDims = ndims(trajectory(acqData))
  if encodingDims==3 && numSlices(acqData)>1
    @error "reconstruction of multiple 3d-encoded volumina is not yet supported"
  end

  if !haskey(recoParams, :reconSize)
    recoParams[:reconSize] = encodingSize(acqData)
  end  

  # load reconstruction parameters
  recoParams = merge(defaultRecoParams(), recoParams)

  # iterative reco
  if recoParams[:reco] == "direct"
    reconSize, weights, cmap = setupDirectReco(acqData, recoParams)
    return reconstruction_direct(acqData, reconSize[1:encodingDims], weights, cmap)
  else
    reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = setupIterativeReco(acqData, recoParams)
    if recoParams[:reco] == "standard"
        return reconstruction_simple(acqData, reconSize[1:encodingDims], reg, sparseTrafo, weights, solvername, normalize, encOps, recoParams)
    elseif recoParams[:reco] == "multiEcho"
        return reconstruction_multiEcho(acqData, reconSize[1:encodingDims], reg, sparseTrafo, weights, solvername, normalize, encOps, recoParams)
    elseif recoParams[:reco] == "multiCoil"
        return reconstruction_multiCoil(acqData, reconSize[1:encodingDims], reg, sparseTrafo, weights, L_inv, solvername, senseMaps, normalize, encOps, recoParams)
    elseif recoParams[:reco] == "multiCoilMultiEcho"
        return reconstruction_multiCoilMultiEcho(acqData, reconSize[1:encodingDims], reg, sparseTrafo, weights, solvername, senseMaps, normalize, encOps, recoParams)
    else
        @error "reco modell $(recoParams[:reco]) not found"
    end
  end
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
