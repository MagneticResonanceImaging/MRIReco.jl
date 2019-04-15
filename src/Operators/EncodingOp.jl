export EncodingOp, lrEncodingOp

#########################
# simple Fourier Encoding
#########################
function encodingOps2d_simple(acqData::AcquisitionData, params::Dict;slice=1)
  tr = [trajectory(acqData,i) for i=1:acqData.numEchoes]
  return [fourierEncodingOp2d(tr[i], params, slice=slice) for i=1:acqData.numEchoes]
end

function encodingOps3d_simple(acqData::AcquisitionData, params::Dict)
  tr = [trajectory(acqData,i) for i=1:acqData.numEchoes]
  return [fourierEncodingOp3d(tr[i], params) for i=1:acqData.numEchoes]
end

##########################################
# fourier encoding with coil sensitivities
##########################################
function encodingOps2d_parallel(acqData::AcquisitionData, params::Dict;slice=1)
  # fourier operators
  ft = encodingOps2d_simple(acqData, params, slice=slice)
  S = SensitivityOp(reshape(params[:senseMaps][:,:,slice,:],:,acqData.numCoils),1)
  return [ diagOp( [ft[i] for k=1:acqData.numCoils]... )*S for i=1:acqData.numEchoes]
end

function encodingOps3d_parallel(acqData::AcquisitionData, params::Dict)
  # fourier operators
  ft = encodingOps3d_simple(acqData, params)
  S = SensitivityOp(reshape(params[:senseMaps],:,acqData.numCoils),1)
  return [ diagOp( [ft[i] for k=1:acqData.numCoils]... )*S for i=1:acqData.numEchoes]
end

#############################
# multi echo fourier encoding
#############################
function encodingOp2d_multiEcho(acqData::AcquisitionData, params::Dict;slice=1)
  # fourier operators
  ft = encodingOps2d_simple(acqData, params, slice=slice)
  return diagOp(ft)
end

function encodingOp3d_multiEcho(acqData::AcquisitionData, params::Dict)
  # fourier operators
  ft = encodingOps3d_simple(acqData, params)
  return diagOp(ft)
end

#####################################################
# multi echo fourier encoding with coil sensitivities
#####################################################
function encodingOp2d_multiEcho_parallel(acqData::AcquisitionData, params::Dict;slice=1)
  # fourier operators
  ft = encodingOps2d_simple(acqData, params, slice=slice)
  S = SensitivityOp(reshape(params[:senseMaps][:,:,slice,:],:,acqData.numCoils),1)
  return diagOp( repeat(ft, acqData.numCoils)... )*S
end

function encodingOp3d_multiEcho_parallel(acqData::AcquisitionData, params::Dict)
  # fourier operators
  ft = encodingOps3d_simple(acqData, params)
  S = SensitivityOp(reshape(params[:senseMaps],:,acqData.numCoils),1)
  return diagOp( repeat(ft, acqData.numCoils)... )*S
end

###################################
# Encoding with low rank projection
###################################
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

  #  TODO: sampling operator in case of undersampled cartesian acquisitions
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
function fourierEncodingOp2d(tr, params; slice=0)
  shape = params[:shape]
  opName="fast"
  if get(params,:fft,false)
    opName="fft"
  elseif get(params,:explicit,false)
    opName="explicit"
  end

  return fourierEncodingOp2d(shape,tr,opName;slice=slice,params...)
end

"""
return 2d Fourier encoding operator (either Explicit, FFT or NFFT)
  opname : "explicit", "fft" or "fast"
  slice : slice to which the operator will be applied
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp2d(shape::NTuple{2,Int64}, tr::Trajectory, opName::String;
          slice::Int64=1, correctionMap=[], echoImage::Bool=true,
          method::String="nfft", alpha::Float64=1.75, m::Float64=4.0, K::Int64=20, kargs...)

  if opName=="explicit"
    @debug "ExplicitOp"
    return ExplicitOp(shape,tr,correctionMap[:,:,slice], echoImage=echoImage)
  elseif opName=="fft"
    @debug "FFTOp"
    return FFTOp(ComplexF64, shape)
  elseif opName=="fast"
    @debug "NFFT-based Op"
    if isempty(correctionMap) || correctionMap==zeros(ComplexF64,size(correctionMap))
      return NFFTOp(shape, tr)
    else
      return FieldmapNFFTOp(shape, tr, correctionMap[:,:,slice],
                            echoImage=echoImage, alpha=alpha, m=m, K=K)
    end
  else
    error("opName $(opName) is not known")
  end

  return NFFTOp(shape, tr)
end

function fourierEncodingOp3d(tr, params)
  shape = params[:shape]
  opName="fast"
  if get(params,:fft,false)
    opName=="fft"
  elseif get(params,:explicit,false)
    opName=="explicit"
  end

  return fourierEncodingOp3d(shape,tr,opName;params...)
end

"""
return 3d Fourier encoding operator (either Explicit, FFT or NFFT)
  opname : "explicit", "fft" or "fast"
  echoImage : calculate signal evolution relative to the echo time
"""
function fourierEncodingOp3d(shape::NTuple{3,Int64}, tr::Trajectory, opName::String;
          correctionMap=[], echoImage::Bool=true,
          method::String="nfft", alpha::Float64=1.75, m::Float64=4.0, K::Int64=20, kargs...)

  if opName=="explicit"
    return ExplicitOp(shape,tr,correctionMap,echoImage=echoImage)
  elseif opName=="fft"
    return FFTOp(ComplexF64, shape)
  elseif opName=="fast"
    if isempty(correctionMap) || correctionMap==zeros(ComplexF64,size(correctionMap))
      return NFFTOp(shape, tr)
    else
      return FieldmapNFFTOp(shape, tr, correctionMap,
                            echoImage=echoImage, alpha=alpha, m=m, K=K)
    end
  else
    error("opName $(opName) is not known")
  end

  return NFFTOp(shape, tr)
end
