export EncodingOp, lrEncodingOp, fourierEncodingOp2d, fourierEncodingOp3d,
       encodingOps2d_simple, encodingOps3d_simple, encodingOps2d_parallel,
       encodingOps3d_parallel, encodingOp2d_multiEcho, encodingOp3d_multiEcho,
       encodingOp2d_multiEcho_parallel, encodingOp3d_multiEcho_parallel

"""
    encodingOps2d_simple(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              ; slice=1
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft")

generates an Array of LinearOperators which describe 2d signal encoding of the individual
contrasts in an MRI acquisition (for a given slice).

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
* `slice=1`                             - slice to be encoded/reconstructed
* (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
"""
function encodingOps2d_simple(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              ; slice=1
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft")

  numContr = numContrasts(acqData)
  tr = [trajectory(acqData,i) for i=1:numContr]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp2d(shape, tr[i], "fast", slice=slice, subsampleIdx=idx[i], correctionMap=correctionMap, method=method) for i=1:numContr]
end

"""
    encodingOps3d_simple(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              ; correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft")

generates an Array of LinearOperators which describe 3d signal encoding of the individual
contrasts&coils in an MRI acquisition (for a given slice).
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOps3d_simple(acqData::AcquisitionData, shape::NTuple{3,Int64}
                              ; correctionMap=ComplexF64[]
                              , method::String="nfft")

  numContr = numContrasts(acqData)
  tr = [trajectory(acqData,i) for i=1:numContr]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp3d(shape, tr[i], "fast", subsampleIdx=idx[i], method=method) for i=1:numContr]
end

"""
    encodingOps2d_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              , senseMaps::Array{ComplexF64}
                              ; slice=1
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft")

generates an Array of LinearOperators which describe 2d signal encoding of the individual
contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
* `slice=1`                             - slice to be encoded/reconstructed
* (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
"""
function encodingOps2d_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                , senseMaps::Array{ComplexF64}
                                ; slice=1
                                , correctionMap::Array{ComplexF64}=ComplexF64[]
                                , method::String="nfft")

  numContr, numChan = numContrasts(acqData), numChannels(acqData)
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape, slice=slice, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps[:,:,slice,:],:,numChan),1)
  return [ diagOp( [ft[i] for k=1:numChan]... )*S for i=1:numContr]
end

"""
    encodingOps3d_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              , senseMaps::Array{ComplexF64}
                              ; correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft")

generates an Array of LinearOperators which describe 3d signal encoding of the individual
contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOps3d_parallel(acqData::AcquisitionData, shape::NTuple{3,Int64}
                                , senseMaps::Array{ComplexF64}
                                ; correctionMap::Array{ComplexF64}=ComplexF64[]
                                , method::String="nfft")

  numContr, numChan = numContrasts(acqData), numChannels(acqData)
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps,:,numChan),1)
  return [ diagOp( [ft[i] for k=1:numChan]... )*S for i=1:numContr]
end

"""
    encodingOp2d_multiEcho(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                ; slice::Int64=1
                                , correctionMap::Array{ComplexF64}=ComplexF64[]
                                , method::String="nfft")

generates a LinearOperator which describe combined 2d signal encoding of all
the contrasts in an MRI acquisition (for a given slice).

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
* `slice=1`                             - slice to be encoded/reconstructed
* (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
"""
function encodingOp2d_multiEcho(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                ; slice::Int64=1
                                , correctionMap::Array{ComplexF64}=ComplexF64[]
                                , method::String="nfft")
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape, slice=slice, correctionMap=correctionMap)
  return diagOp(ft...)
end

"""
    encodingOp3d_multiEcho(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                ; correctionMap::Array{ComplexF64}=ComplexF64[])

generates a LinearOperator which describe combined 3d signal encoding of all
the contrasts in an MRI acquisition (for a given slice).
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOp3d_multiEcho(acqData::AcquisitionData, shape::NTuple{3,Int64}; correctionMap::Array{ComplexF64,3}=ComplexF64[])
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape, correctionMap=correctionMap)
  return diagOp(ft...)
end

"""
    encodingOp2d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; slice=1
                                          , correctionMap::Array{ComplexF64}=ComplexF64[]
                                          , method::String="nfft")

generates a LinearOperator which describe combined 3d signal encoding of all
the contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
* `slice=1`                             - slice to be encoded/reconstructed
* (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
"""
function encodingOp2d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; slice=1
                                          , correctionMap::Array{ComplexF64}=ComplexF64[]
                                          , method::String="nfft")

  numChan = numChannels(acqData)
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape, slice=slice, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps[:,:,slice,:],:,numChan),numContrasts(acqData))
  return diagOp( repeat(ft, numChan)... )*S
end

"""
    encodingOp3d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; correctionMap::Array{ComplexF64}=ComplexF64[]
                                          , method::String="nfft")

generates a LinearOperator which describe combined 3d signal encoding of all
the contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities
Arguments are the same as in the 2d case, with the exception that shape is of type `NTuple{3,Int64}`
and the considered slice is not specified.
"""
function encodingOp3d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{3,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; correctionMap::Array{ComplexF64}=ComplexF64[]
                                          , method::String="nfft")

  numChan = numChannels(acqData)
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps,:,numChan),1)
  return diagOp( repeat(ft, numChan)... )*S
end

###################################
# Encoding with low rank projection
###################################
function lrEncodingOp(acqData::AcquisitionData, shape, params::Dict; numContr::Int64=1, parallel::Bool=false,)

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

  return M*Φ*E

end

"""
return 2d Fourier encoding operator (either Explicit or NFFT)
  * `opname` : "explicit" or "fast"
  * `slice` : slice to which the operator will be applied
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp2d(shape::NTuple{2,Int64}, tr::Trajectory, opName::String;
          subsampleIdx::Vector{Int64}=Int64[], slice::Int64=1, correctionMap=[]
          , echoImage::Bool=true, method::String="nfft", kargs...)
  # Fourier transformations
  if opName=="explicit"
    @debug "ExplicitOp"
    ftOp = ExplicitOp(shape,tr,correctionMap[:,:,slice], echoImage=echoImage)
  elseif opName=="fast"
    @debug "NFFT-based Op"
    if !isempty(correctionMap) && correctionMap!=zeros(ComplexF64,size(correctionMap))
      ftOp = FieldmapNFFTOp(shape, tr, correctionMap[:,:,slice], echoImage=echoImage, method=method)
    elseif isCartesian(tr)
      @debug "FFTOp"
      ftOp = sqrt(prod(shape))*FFTOp(ComplexF64, shape)
    else
      ftOp = NFFTOp(shape, tr)
    end
  else
    @error "opName $(opName) is not known"
  end

  # subsampling
  if !isempty(subsampleIdx) && length(subsampleIdx)!=size(tr,2)
    S = SamplingOp(subsampleIdx,shape)
  else
    S = opEye()
  end

  return S*ftOp
end

"""
return 3d Fourier encoding operator (either Explicit, FFT or NFFT)
  opname : "explicit", "fft" or "fast"
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp3d(shape::NTuple{3,Int64}, tr::Trajectory, opName::String
          ; subsampleIdx::Vector{Int64}=Int64[], correctionMap=ComplexF64[], echoImage::Bool=true
          , kargs...)

  if opName=="explicit"
    ftOp = ExplicitOp(shape,tr,correctionMap,echoImage=echoImage)
  elseif opName=="fast"
    if !isempty(correctionMap) && correctionMap!=zeros(ComplexF64,size(correctionMap))
      ftOp = FieldmapNFFTOp(shape, tr, correctionMap, echoImage=echoImage)
    elseif isCartesian(tr)
      ftOp = sqrt(prod(shape))*FFTOp(ComplexF64, shape)
    else
      ftOp = NFFTOp(shape, tr)
    end
  else
    error("opName $(opName) is not known")
  end

  # subsampling
  if !isempty(subsampleIdx) && length(subsampleIdx)!=size(tr,2)
    S = SamplingOp(subsampleIdx,shape)
  else
    S = opEye()
  end

  return S*ftOp
end
