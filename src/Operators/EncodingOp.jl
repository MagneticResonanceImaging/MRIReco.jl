export EncodingOp, lrEncodingOp

#
# Encoding operator for one slice
#
function EncodingOp2d(acqData::AcquisitionData, params::Dict; slice=1,
                           parallel=false, multiEcho=false)
  # Fourier Operators
  ft = Array{AbstractLinearOperator,1}(undef,acqData.numEchoes)

  for i = 1:acqData.numEchoes
    tr = trajectory(acqData,i)
    ft[i] = fourierEncodingOp2d(tr, params, slice=slice)
  end

  # coil sensitivities in case of SENSE-reconstruction
  if parallel && multiEcho
    S = SensitivityOp( reshape(params[:senseMaps][:,:,slice,:],:,acqData.numCoils),
                       acqData.numEchoes)
    E = diagOp( repeat(ft, acqData.numCoils)... )*S
  elseif parallel && !multiEcho
    S = SensitivityOp(reshape(params[:senseMaps][:,:,slice,:],:,acqData.numCoils),1)
    E = [ diagOp( [ft[i] for k=1:acqData.numCoils]... )*S for i=1:acqData.numEchoes]
  elseif !parallel && multiEcho
    E = diagOp(ft...)
  else
    E = ft
  end

  # sampling operator in case of undersampled cartesian acquisitions
  if acqData.subsampleIndices != collect(1:length(acqData.kdata))
    if get(params,:sampling, nothing) == "binary"
      M = SamplingOp(Array{Bool})(acqData.subsampleIndices)
    else
      numNodes = prod(params[:shape])
      M = SamplingOp(acqData.subsampleIndices,(numNodes, acqData.numEchoes, acqData.numCoils))
    end
    return M*E
  end

  return E
end

function EncodingOp3d(acqData::AcquisitionData, params::Dict;
                           parallel=false, multiEcho=false)
  # Fourier Operators
  ft = Array{AbstractLinearOperator,1}(undef,acqData.numEchoes)

  for i = 1:acqData.numEchoes
    tr = trajectory(acqData,i)
    ft[i] = fourierEncodingOp3d(tr, params)
  end

  # coil sensitivities in case of SENSE-reconstruction
  if parallel && multiEcho
    S = SensitivityOp( reshape(params[:senseMaps],:,acqData.numCoils),
                       acqData.numEchoes)
    E = diagOp( repeat(ft, acqData.numCoils)... )*S
  elseif parallel && !multiEcho
    S = SensitivityOp(reshape(params[:senseMaps],:,acqData.numCoils),1)
    E = [ diagOp( [ft[i] for k=1:acqData.numCoils]... )*S for i=1:acqData.numEchoes]
  elseif !parallel && multiEcho
    E = diagOp(ft...)
  else
    E = ft
  end

  # sampling operator in case of undersampled cartesian acquisitions
  if acqData.subsampleIndices != collect(1:length(acqData.kdata))
    if get(params,:sampling, nothing) == "binary"
      M = SamplingOp(Array{Bool})(acqData.subsampleIndices)
    else
      numNodes = prod(params[:shape])
      M = SamplingOp(acqData.subsampleIndices,(numNodes, acqData.numEchoes, acqData.numCoils))
    end
    return M*E
  end

  return E
end

function lrEncodingOp(acqData::AcquisitionData,params::Dict; numEchoes::Int64=1, parallel::Bool=false,)

  # low rank operator
  N=prod(params[:shape])
  subspace = get( params, :phi, Matrix{Float64}(I,N,N) )
  K = size(subspace,2)
  Φ = MapSliceOp(subspace[:,:,1],2,(N, K, acqData.numCoils), (N, numEchoes, acqData.numCoils))

  # Fourier Operator
  tr = trajectory(acqData,1)
  params[:fft] = true
  ft = fourierEncodingOp2d(tr, params)

  # coil sensitivities in case of SENSE-reconstruction
  if parallel
    S = SensitivityOp(params[:senseMaps][:,:],K)
    E = diagOp( [ft for i=1:acqData.numCoils*K]... )*S
  else
    E = diagOp( [ft for i=1:K]... )
  end

  # sampling operator in case of undersampled cartesian acquisitions
  if get(params,:sampling, nothing) == "binary"
    M = SamplingOp(Array{Bool}(acqData.subsampleIndices))
  else
    M = SamplingOp(acqData.subsampleIndices,(N, numEchoes, acqData.numCoils))
  end

  return M*Φ*E

end


#
# return Fourier encoding operator with(out) correction when cmap is (not) specified
#
# function fourierEncodingOp(tr, params; slice=1)
#   densityWeighting = get(params,:densityWeighting,true)
#   fft = get(params,:fft,false)
#   if fft
#     FFTOp(ComplexF64, params[:shape])
#   elseif !haskey(params,:correctionMap)
#     return NFFTOp(params[:shape], tr, symmetrize=densityWeighting)
#   else
#     echoImage = get(params, :echoImage, true)
#     method = get(params, :method, "nfft")
#     alpha = get(params, :alpha, 1.75)
#     m = get(params, :m, 1.75)
#     K = get(params, :K, 20)
#     return FieldmapNFFTOp( params[:shape],  tr, params[:correctionMap][:,:,slice],
#                            symmetrize=densityWeighting, echoImage=echoImage, K=K,
#                            alpha=alpha, m=m, method=method )
#   end
# end

function fourierEncodingOp2d(tr, params; slice=0)
  shape = params[:shape]
  opName="fast"
  if get(params,:fft,false)
    opName="fft"
  elseif get(params,:explicit,false)
    opName="explicit"
  end
  symmetrize = get(params, :densityWeighting, true)

  return fourierEncodingOp2d(shape,tr,opName;slice=slice, symmetrize=symmetrize,params...)
end

"""
return 2d Fourier encoding operator (either Explicit, FFT or NFFT)
  opname : "explicit", "fft" or "fast"
  slice : slice to which the operator will be applied
  symmetrize : symmetrizes the operator
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp2d(shape::NTuple{2,Int64}, tr::Trajectory, opName::String;
          slice::Int64=1, correctionMap=[], symmetrize::Bool=true, echoImage::Bool=true,
          method::String="nfft", alpha::Float64=1.75, m::Float64=4.0, K::Int64=20, kargs...)

  if opName=="explicit"
    @debug "ExplicitOp"
    return ExplicitOp(shape,tr,correctionMap[:,:,slice], symmetrize=symmetrize, echoImage=echoImage)
  elseif opName=="fft"
    @debug "FFTOp"
    return FFTOp(ComplexF64, shape)
  elseif opName=="fast"
    @debug "NFFT-based Op"
    if isempty(correctionMap) || correctionMap==zeros(ComplexF64,size(correctionMap))
      return NFFTOp(shape, tr, symmetrize=symmetrize)
    else
      return FieldmapNFFTOp(shape, tr, correctionMap[:,:,slice], symmetrize=symmetrize,
                            echoImage=echoImage, alpha=alpha, m=m, K=K)
    end
  else
    error("opName $(opName) is not known")
  end

  return NFFTOp(shape, tr, symmetrize=symmetrize)
end

function fourierEncodingOp3d(tr, params)
  shape = params[:shape]
  opName="fast"
  if get(params,:fft,false)
    opName=="fft"
  elseif get(params,:explicit,false)
    opName=="explicit"
  end
  symmetrize = get(params, :densityWeighting, true)

  return fourierEncodingOp3d(shape,tr,opName; symmetrize=symmetrize,params...)
end

"""
return 3d Fourier encoding operator (either Explicit, FFT or NFFT)
  opname : "explicit", "fft" or "fast"
  symmetrize : symmetrizes the operator
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp3d(shape::NTuple{3,Int64}, tr::Trajectory, opName::String;
          correctionMap=[], symmetrize::Bool=true, echoImage::Bool=true,
          method::String="nfft", alpha::Float64=1.75, m::Float64=4.0, K::Int64=20, kargs...)

  if opName=="explicit"
    return ExplicitOp(shape,tr,correctionMap, symmetrize=symmetrize, echoImage=echoImage)
  elseif opName=="fft"
    return FFTOp(ComplexF64, shape)
  elseif opName=="fast"
    if isempty(correctionMap) || correctionMap==zeros(ComplexF64,size(correctionMap))
      return NFFTOp(shape, tr, symmetrize=symmetrize)
    else
      return FieldmapNFFTOp(shape, tr, correctionMap, symmetrize=symmetrize,
                            echoImage=echoImage, alpha=alpha, m=m, K=K)
    end
  else
    error("opName $(opName) is not known")
  end

  return NFFTOp(shape, tr, symmetrize=symmetrize)
end
