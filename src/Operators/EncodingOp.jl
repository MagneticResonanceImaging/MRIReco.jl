export EncodingOp, lrEncodingOp, fourierEncodingOp2d, fourierEncodingOp3d,
       encodingOps2d_simple, encodingOps3d_simple, encodingOps2d_parallel,
       encodingOps3d_parallel, encodingOp2d_multiEcho, encodingOp3d_multiEcho,
       encodingOp2d_multiEcho_parallel, encodingOp3d_multiEcho_parallel

"""
    encodingOps2d_simple(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              ; kargs...)

generates an Array of LinearOperators which describe 2d signal encoding of the individual
contrasts in an MRI acquisition (for a given slice).

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
"""
function encodingOps2d_simple(acqData::AcquisitionData, shape::NTuple{2,Int64}; kargs...)
  numContr = numContrasts(acqData)
  tr = [trajectory(acqData,i) for i=1:numContr]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp2d(shape, tr[i], "fast", 
                 subsampleIdx=idx[i]; kargs...) for i=1:numContr]
end

"""
    encodingOps3d_simple(acqData::AcquisitionData, shape::NTuple{2,Int64}; kargs...)

generates an Array of LinearOperators which describe 3d signal encoding of the individual
contrasts&coils in an MRI acquisition (for a given slice).
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOps3d_simple(acqData::AcquisitionData, shape::NTuple{3,Int64}; kargs...)
  numContr = numContrasts(acqData)
  tr = [trajectory(acqData,i) for i=1:numContr]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp3d(shape, tr[i], "fast", subsampleIdx=idx[i]; kargs...) for i=1:numContr]
end

"""
    encodingOps2d_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              , senseMaps::Array{ComplexF64}
                              ; kargs...)

generates an Array of LinearOperators which describe 2d signal encoding of the individual
contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
"""
function encodingOps2d_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                , senseMaps::Array{ComplexF64}
                                ; slice=1, kargs...)

  numContr, numChan = numContrasts(acqData), numChannels(acqData)
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape; slice=slice, kargs...)
  S = SensitivityOp(reshape(senseMaps[:,:,slice,:],:,numChan),1)
  Op = [ diagOp(ft[i], numChan) ∘ S for i=1:numContr]

  return Op
end

"""
    encodingOps3d_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              , senseMaps::Array{ComplexF64}
                              ; kargs...)

generates an Array of LinearOperators which describe 3d signal encoding of the individual
contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOps3d_parallel(acqData::AcquisitionData, shape::NTuple{3,Int64}
                                , senseMaps::Array{ComplexF64}
                                ; kargs...)

  numContr, numChan = numContrasts(acqData), numChannels(acqData)
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape; kargs...)
  S = SensitivityOp(reshape(senseMaps,:,numChan),1)
  return [ diagOp(ft[i], numChan) ∘ S for i=1:numContr]
end

"""
    encodingOp2d_multiEcho(acqData::AcquisitionData, shape::NTuple{2,Int64}; kargs...)

generates a LinearOperator which describe combined 2d signal encoding of all
the contrasts in an MRI acquisition (for a given slice).

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
"""
function encodingOp2d_multiEcho(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                ; kargs...)
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape; kargs...)
  return diagOp(ft...)
end

"""
    encodingOp3d_multiEcho(acqData::AcquisitionData, shape::NTuple{2,Int64}; kargs...)

generates a LinearOperator which describe combined 3d signal encoding of all
the contrasts in an MRI acquisition (for a given slice).
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOp3d_multiEcho(acqData::AcquisitionData, shape::NTuple{3,Int64}; kargs...)
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape; kargs...)
  return diagOp(ft...)
end

"""
    encodingOp2d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; kargs...)

generates a LinearOperator which describe combined 3d signal encoding of all
the contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
"""
function encodingOp2d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; slice::Int64=1, kargs...)

  numChan = numChannels(acqData)
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape; kargs...)
  S = SensitivityOp(reshape(senseMaps[:,:,slice,:],:,numChan),numContrasts(acqData))
  return diagOp(ft, numChan) ∘ S
end

"""
    encodingOp3d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; kargs...)

generates a LinearOperator which describe combined 3d signal encoding of all
the contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOp3d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{3,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; kargs...)

  numChan = numChannels(acqData)
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape; kargs...)
  S = SensitivityOp(reshape(senseMaps,:,numChan),1)
  return diagOp(ft, numChan) ∘ S
end

###################################
# Encoding with low rank projection
###################################
function lrEncodingOp(acqData::AcquisitionData, shape, params::Dict; numContr::Int64=1, parallel::Bool=false)

  numChan = numChannels(acqData)
  # low rank operator
  N=prod(shape)
  subspace = get( params, :phi, Matrix{Float64}(I,N,N) )
  K = size(subspace,2)
  Φ = MapSliceOp(subspace[:,:,1],2,(N, K, numChan), (N, numContr, numChan))

  # Fourier Operator
  tr = trajectory(acqData,1)
  ft = fourierEncodingOp2d(shape, tr, "fast"; params...)

  # coil sensitivities in case of SENSE-reconstruction
  if parallel
    S = SensitivityOp(params[:senseMaps][:,:],K)
    E = diagOp( [ft for i=1:numChan*K]... )*S
  else
    E = diagOp( [ft for i=1:K]... )
  end

  #  TODO: sampling operator in case of undersampled cartesian acquisitions
  # vectorize sampling pattern
  if get(params,:sampling, nothing) == "binary"
    M = SamplingOp(Array{Bool}(acqData.subsampleIndices))
  else
    # sampling idx for all contrasts but one coil
    subIdx = hcat(acqData.subsampleIndices...)
    M = SamplingOp( hcat([subIdx for c=1:numChan]...), (N, numContr, numChan))
  end

  return M ∘ (Φ ∘ E)

end

"""
return 2d Fourier encoding operator (either Explicit or NFFT)
  * `opname` : "explicit" or "fast"
  * `slice` : slice to which the operator will be applied
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp2d(shape::NTuple{2,Int64}, tr::Trajectory, opName::String;
          subsampleIdx::Vector{Int64}=Int64[], slice::Int64=1, correctionMap=[],
          echoImage::Bool=true,  kargs...)
  # Fourier transformations
  if opName=="explicit"
    @debug "ExplicitOp"
    ftOp = ExplicitOp(shape,tr,correctionMap[:,:,slice], echoImage=echoImage)
  elseif opName=="fast"
    @debug "NFFT-based Op"
    if !isempty(correctionMap) && correctionMap!=zeros(ComplexF64,size(correctionMap))
      ftOp = FieldmapNFFTOp(shape, tr, correctionMap[:,:,slice], echoImage=echoImage; kargs...)
    elseif isCartesian(tr)
      @debug "FFTOp"
      ftOp = FFTOp(ComplexF64, shape; unitary=false)
    else
      ftOp = NFFTOp(shape, tr; kargs...)
    end
  else
    @error "opName $(opName) is not known"
  end

  # subsampling
  if !isempty(subsampleIdx) && (subsampleIdx != collect(1:size(tr,2))) && isCartesian(tr)
    S = SamplingOp(subsampleIdx,(tr.numSamplingPerProfile,tr.numProfiles))
    return S ∘ ftOp
  else
    return ftOp
  end
end

"""
return 3d Fourier encoding operator (either Explicit, FFT or NFFT)
  opname : "explicit", "fft" or "fast"
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp3d(shape::NTuple{3,Int64}, tr::Trajectory, opName::String;
          subsampleIdx::Vector{Int64}=Int64[], correctionMap=ComplexF64[], echoImage::Bool=true,
          kargs...)

  if opName=="explicit"
    ftOp = ExplicitOp(shape,tr,correctionMap,echoImage=echoImage)
  elseif opName=="fast"
    if !isempty(correctionMap) && correctionMap!=zeros(ComplexF64,size(correctionMap))
      ftOp = FieldmapNFFTOp(shape, tr, correctionMap, echoImage=echoImage; kargs...)
    elseif isCartesian(tr)
      #ftOp = sqrt(prod(shape))*FFTOp(ComplexF64, shape)
      ftOp = FFTOp(ComplexF64, shape; unitary=false)
    else
      ftOp = NFFTOp(shape, tr; kargs...)
    end
  else
    error("opName $(opName) is not known")
  end

  # subsampling
  if !isempty(subsampleIdx) && (subsampleIdx != collect(1:size(tr,2))) && isCartesian(tr)
    S = SamplingOp(subsampleIdx,shape)
    return S ∘ ftOp
  else
    return ftOp
  end
end
