export EncodingOp, lrEncodingOp, fourierEncodingOp, encodingOps_simple, 
       encodingOps_parallel, encodingOp_multiEcho, encodingOp_multiEcho_parallel

"""
    encodingOps_simple(acqData::AcquisitionData, shape::NTuple{D,Int64}
                              ; kargs...) where D

generates an Array of LinearOperators which describe signal encoding of the individual
contrasts in an MRI acquisition (for a given slice).

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{D,Int64}`              - size of image to be encoded/reconstructed
"""
function encodingOps_simple(acqData::AcquisitionData{T}, shape::NTuple{D,Int64}; kargs...) where {T,D}
  numContr = numContrasts(acqData)
  tr = [trajectory(acqData,i) for i=1:numContr]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp(shape, tr[i], "fast", 
                 subsampleIdx=idx[i]; kargs...) for i=1:numContr]
end

"""
    encodingOps_parallel(acqData::AcquisitionData{T}, shape::NTuple{D,Int64}
                              , senseMaps::Array{Complex{T}}
                              ; kargs...) where {T,D}

generates an Array of LinearOperators which describe signal encoding of the individual
contrasts in an MRI acquisition. The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData{T}`            - AcquisitionData object
* `shape::NTuple{D,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::Array{Complex{T}}`        - coil sensitivities
"""
function encodingOps_parallel(acqData::AcquisitionData{T}, shape::NTuple{D,Int64}
                                , senseMaps::Array{Complex{T},4}
                                ; slice=1, kargs...) where {T,D}

  smaps = ( D==2 ? senseMaps[:,:,slice,:] : senseMaps )

  numContr, numChan = numContrasts(acqData), numChannels(acqData)
  # fourier operators
  ft = encodingOps_simple(acqData, shape; slice=slice, kargs...)
  S = SensitivityOp(reshape(smaps,:,numChan),1)
  Op = [ diagOp(ft[i], numChan) ∘ S for i=1:numContr]

  return Op
end

"""
    encodingOp_multiEcho(acqData::AcquisitionData, shape::NTuple{D,Int64}; kargs...) where D

generates a LinearOperator which describe combined signal encoding of all
the contrasts in an MRI acquisition (for a given slice).

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{D,Int64}`              - size of image to be encoded/reconstructed
"""
function encodingOp_multiEcho(acqData::AcquisitionData, shape::NTuple{D,Int64}
                                ; kargs...) where D
  # fourier operators
  ft = encodingOps_simple(acqData, shape; kargs...)
  return diagOp(ft...)
end

"""
    encodingOp_multiEcho_parallel(acqData::AcquisitionData{T}, shape::NTuple{D,Int64}
                                          , senseMaps::Array{Complex{T}}
                                          ; kargs...) where {T,D}

generates a LinearOperator which describe combined signal encoding of all
the contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData{T}`            - AcquisitionData object
* `shape::NTuple{2,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::Array{Complex{T}}`        - coil sensitivities
"""
function encodingOp_multiEcho_parallel(acqData::AcquisitionData{T}, shape::NTuple{D,Int64}
                                          , senseMaps::Array{Complex{T}}
                                          ; slice::Int64=1, kargs...) where {T,D}

  smaps = ( D==2 ? senseMaps[:,:,slice,:] : senseMaps )

  numChan = numChannels(acqData)
  # fourier operators
  ft = encodingOps_simple(acqData, shape; kargs...)
  S = SensitivityOp(reshape(smaps,:,numChan),numContrasts(acqData)) 
  return diagOp(repeat(ft, outer=numChan)...) ∘ S
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
return Fourier encoding operator (either Explicit or NFFT)
  * `opname` : "explicit" or "fast"
  * `slice` : slice to which the operator will be applied
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp(shape::NTuple{D,Int64}, tr::Trajectory{T}, opName::String;
          subsampleIdx::Vector{Int64}=Int64[], slice::Int64=1, correctionMap::Array{Complex{T}}=Complex{T}[],
          echoImage::Bool=true,  kargs...) where {T,D}

  # extract proper portion of correctionMap
  if !isempty(correctionMap)
    cmap = ( D==2 ? correctionMap[:,:,slice] : correctionMap )
  end
  # Fourier transformations
  if opName=="explicit"
    @debug "ExplicitOp"
    ftOp = ExplicitOp(shape, tr, cmap, echoImage=echoImage)
  elseif opName=="fast"
    @debug "NFFT-based Op"
    if !isempty(correctionMap) && correctionMap!=zeros(Complex{T},size(correctionMap))
      ftOp = FieldmapNFFTOp(shape, tr, cmap, echoImage=echoImage; kargs...)
    elseif isCartesian(tr)
      @debug "FFTOp"
      ftOp = FFTOp(Complex{T}, shape; unitary=false)
    else
      ftOp = NFFTOp(shape, tr; kargs...)
    end
  else
    @error "opName $(opName) is not known"
  end

  # subsampling
  if !isempty(subsampleIdx) && (subsampleIdx != collect(1:size(tr,2))) && isCartesian(tr)
    if D==2
      S = SamplingOp(subsampleIdx,(tr.numSamplingPerProfile, tr.numProfiles), Complex{T})
    else
      S = SamplingOp(subsampleIdx,(tr.numSamplingPerProfile, tr.numProfiles, tr.numSlices), Complex{T})
    end
    return S ∘ ftOp
  else
    return ftOp
  end
end
