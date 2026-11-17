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
function encodingOps_simple(acqData::AcquisitionData{T,D}, shape::NTuple{D,Int64}; kargs...) where {T,D}
  numContr = numContrasts(acqData)
  tr = [trajectory(acqData,i) for i=1:numContr]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp(shape, tr[i], "fast",
                 subsampleIdx=idx[i]; kargs...) for i=1:numContr]
end

"""
    encodingOps_parallel(acqData::AcquisitionData{T,D}, shape::NTuple{D,Int64}
                              , senseMaps::AbstractArray{Complex{T},4}
                              ; kargs...) where {T,D}

generates an Array of LinearOperators which describe signal encoding of the individual
contrasts in an MRI acquisition. The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData{T,D}`       - AcquisitionData object
* `shape::NTuple{D,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::AbstractArray{Complex{T},4}` - coil sensitivities with the coils in the fourth dimension
"""
function encodingOps_parallel(acqData::AcquisitionData{T,D}, shape::NTuple{D,Int64}
                                , senseMaps::AbstractArray{Complex{T},4}
                                ; slice=1, S = Vector{Complex{T}}, copyOpsFn = copy, kargs...) where {T,D}

  smaps = ( D==2 ? senseMaps[:,:,slice,:] : senseMaps )

  numContr, numChan = numContrasts(acqData), numChannels(acqData)
  # fourier operators
  ft = encodingOps_simple(acqData, shape; slice=slice, S = S, kargs...)
  SOp = SensitivityOp(reshape(smaps,:,numChan),1)
  Op = [ DiagOp(ft[i], numChan; copyOpsFn = copyOpsFn) ∘ SOp for i=1:numContr]

  return Op
end

"""
    encodingOp_multiEcho(acqData::AcquisitionData{T,D}, shape::NTuple{D,Int64}; kargs...) where {T,D}

generates a LinearOperator which describe combined signal encoding of all
the contrasts in an MRI acquisition (for a given slice).

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `shape::NTuple{D,Int64}`              - size of image to be encoded/reconstructed
"""
function encodingOp_multiEcho(acqData::AcquisitionData{T,D}, shape::NTuple{D,Int64}
                                ; kargs...) where {T,D}
  # fourier operators
  ft = encodingOps_simple(acqData, shape; kargs...)
  return DiagOp(ft...)
end

"""
    encodingOp_multiEcho_parallel(acqData::AcquisitionData{T,D}, shape::NTuple{D,Int64}
                                          , senseMaps::AbstractArray{Complex{T}}
                                          ; kargs...) where {T,D}

generates a LinearOperator which describe combined signal encoding of all
the contrasts in an MRI acquisition (for a given slice). The different coils are taken into account
in terms of their sensitivities

# Arguments
* `acqData::AcquisitionData{T,D}`       - AcquisitionData object
* `shape::NTuple{D,Int64}`              - size of image to be encoded/reconstructed
* `senseMaps::AbstractArray{Complex{T}}` - coil sensitivities with the coils in the fourth dimension
"""
function encodingOp_multiEcho_parallel(acqData::AcquisitionData{T,D}, shape::NTuple{D,Int64}
                                          , senseMaps::AbstractArray{Complex{T}}
                                          ; slice::Int64=1, copyOpsFn = copy, kargs...) where {T,D}

  smaps = ( D==2 ? senseMaps[:,:,slice,:] : senseMaps )

  numChan = numChannels(acqData)
  # fourier operators
  ft = encodingOps_simple(acqData, shape; kargs...)
  SOp = SensitivityOp(reshape(smaps,:,numChan),numContrasts(acqData))
  ops2 = [copyOpsFn(ft[n]) for j=1:numChan,n=eachindex(ft)]
  return DiagOp(ops2...) ∘ SOp
end

###################################
# Encoding with low rank projection
###################################
function lrEncodingOp(acqData::AcquisitionData, shape, params::Dict; numContr::Int64=1, parallel::Bool=false)

  numChan = numChannels(acqData)
  # low rank operator
  N = prod(shape)
  subspace = get( params, :phi, Matrix{Float64}(I,N,N) )
  K = size(subspace,2)
  Φ = MapSliceOp(subspace[:,:,1], 2, (N, K, numChan), (N, numContr, numChan))

  # Fourier Operator
  tr = trajectory(acqData,1)
  ft = fourierEncodingOp(shape, tr, "fast"; params...)

  # coil sensitivities in case of SENSE-reconstruction
  if parallel
    S = SensitivityOp(params[:senseMaps][:,:],K)
    E = DiagOp( [ft for i=1:numChan*K]... )*S
  else
    E = DiagOp( [ft for i=1:K]... )
  end

  #  TODO: sampling operator in case of undersampled cartesian acquisitions
  # vectorize sampling pattern
  if get(params,:sampling, nothing) == "binary"
    M = SamplingOp(ComplexF64, pattern=Array{Bool}(acqData.subsampleIndices))
  else
    # sampling idx for all contrasts but one coil
    subIdx = hcat(acqData.subsampleIndices...)
    M = SamplingOp(ComplexF64, pattern=hcat([subIdx for c=1:numChan]...), shape=(N, numContr, numChan))
  end

  return M ∘ (Φ ∘ E)
end

"""
    fourierEncodingOp(shape::NTuple{D,Int64}, tr::Trajectory{T}, opName::String;
                      subsampleIdx::Vector{Int64}=Int64[], slice::Int64=1,
                      correctionMap::Array{Complex{T}}=Complex{T}[], echoImage::Bool=true,
                      S = Vector{Complex{T}}, kargs...) where {T,D}

returns the Fourier encoding operator of the trajectory `tr` for images of size `shape`.

# Arguments
* `shape::NTuple{D,Int64}`          - size of image to encode/reconstruct
* `tr::Trajectory`                  - trajectory with the k-space nodes to sample
* `opName::String`                  - `"explicit"` for an [`ExplicitOp`](@ref), or `"fast"` for a
                                      [`FieldmapNFFTOp`](@ref) if the correction map is not zero (two-dimensional
                                      images only), otherwise an `NFFTOp` for non-Cartesian trajectories and, for
                                      Cartesian trajectories with no more samples, profiles or slices than `shape`,
                                      an `FFTOp`, with a `SamplingOp` if they have fewer; other names throw an error
* (`subsampleIdx::Vector{Int64}`)   - for Cartesian trajectories, the indices of the sampled nodes; unless the
                                      vector is empty or equals `1:N` for the `N` nodes of `tr`, a `SamplingOp`
                                      selects them after the transform
* (`slice::Int64=1`)                - for two-dimensional images, the slice of the correction map, along its third
                                      dimension, that the operator uses
* (`correctionMap`)                 - complex correction map with the relaxation rate R₂* in 1/s as its real part and
                                      the off-resonance frequency in rad/s as its imaginary part, in the precision
                                      of the trajectory
* (`echoImage::Bool=true`)          - passed to `ExplicitOp` and `FieldmapNFFTOp`
* (`S = Vector{Complex{T}}`)        - storage vector type, for example `CuVector{ComplexF32}` on a GPU
* (`kargs`)                         - options of the approximation for `FieldmapNFFTOp`, and NFFT options such as
                                      `oversamplingFactor` and `kernelSize` for `NFFTOp`, which also takes an
                                      integer `m` as its kernel size; `ExplicitOp` and `FFTOp` ignore them
"""
function fourierEncodingOp(shape::NTuple{D,Int64}, tr::Trajectory{T}, opName::String;
          subsampleIdx::Vector{Int64}=Int64[], slice::Int64=1, correctionMap::Array{Complex{T}}=Complex{T}[],
          echoImage::Bool=true, S = Vector{Complex{T}}, copyOpsFn = identity, kargs...) where {T,D}

  # copyOpsFn is unused in operators, but can break dispatch if passed as keyword argument to operator

  # extract proper portion of correctionMap
  if !isempty(correctionMap)
    cmap = ( D==2 ? correctionMap[:,:,slice] : correctionMap )
  end
  # Fourier transformations
  if opName=="explicit"
    @debug "ExplicitOp"
    if isempty(correctionMap)
      cmap = zeros(Complex{T}, shape)
    end
    ftOp = ExplicitOp(shape, tr, cmap, echoImage=echoImage, S = S)
  elseif opName=="fast"
    @debug "NFFT-based Op"
    if !isempty(correctionMap) && correctionMap!=zeros(Complex{T},size(correctionMap))
      ftOp = FieldmapNFFTOp(shape, tr, cmap, echoImage=echoImage; S = S, fftParams(S)..., kargs...)
    elseif isCartesian(tr)
      @debug "FFTOp"
      if !MRIBase.isUndersampledCartTrajectory(shape,tr)
        ftOp = FFTOp(Complex{T}; shape, unitary=false, S = S, fftParams(S)...)
      else
        idx = MRIBase.cartesianSubsamplingIdx(shape,tr)
        ftOp = SamplingOp(Complex{T}; pattern=idx, shape, S = S) ∘ FFTOp(Complex{T}; shape, unitary=false, S = S, fftParams(S)...)
      end
    else
      # options of the field map approximation are only known to FieldmapNFFTOp
      nfftKargs = filter(kv -> !(first(kv) in (:method, :alpha, :K, :K_tol, :step, :scheduler)), kargs)
      ftOp = NFFTOp(Complex{T}; nodes = kspaceNodes(tr), shape, S = S, nfftParams(S)..., nfftKargs...)
    end
  else
    @error "opName $(opName) is not known"
  end

  # subsampling
  if !isempty(subsampleIdx) && (subsampleIdx != collect(1:size(tr,2))) && isCartesian(tr)
    β = (D==2) ? (tr.numSamplingPerProfile, tr.numProfiles) : (tr.numSamplingPerProfile, tr.numProfiles, tr.numSlices)
    S = SamplingOp(Complex{T}; pattern = subsampleIdx, shape=β, S = S)
    return S ∘ ftOp
  else
    return ftOp
  end
end
