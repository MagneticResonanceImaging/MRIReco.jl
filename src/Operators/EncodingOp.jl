export EncodingOp, lrEncodingOp, fourierEncodingOp2d, fourierEncodingOp3d

#########################
# simple Fourier Encoding
#########################
function encodingOps2d_simple(acqData::AcquisitionData, shape::NTuple{2,Int64}
                              ; slice=1
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft")

  tr = [trajectory(acqData,i) for i=1:acqData.numEchoes]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp2d(shape, tr[i], "fast", slice=slice, subsampleIdx=idx[i], correctionMap=correctionMap, method=method) for i=1:acqData.numEchoes]
end

function encodingOps3d_simple(acqData::AcquisitionData, shape::NTuple{3,Int64}
                              ; correctionMap=ComplexF64[]
                              , method::String="nfft")

  tr = [trajectory(acqData,i) for i=1:acqData.numEchoes]
  idx = acqData.subsampleIndices
  return [fourierEncodingOp3d(shape, tr[i], "fast", subsampleIdx=idx[i], method=method) for i=1:acqData.numEchoes]
end

##########################################
# fourier encoding with coil sensitivities
##########################################
function encodingOps2d_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                , senseMaps::Array{ComplexF64}
                                ; slice=1
                                , correctionMap::Array{ComplexF64}=ComplexF64[]
                                , method::String="nfft")
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape, slice=slice, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps[:,:,slice,:],:,acqData.numCoils),1)
  return [ diagOp( [ft[i] for k=1:acqData.numCoils]... )*S for i=1:acqData.numEchoes]
end

function encodingOps3d_parallel(acqData::AcquisitionData, shape::NTuple{3,Int64}
                                , senseMaps::Array{ComplexF64}
                                ; correctionMap::Array{ComplexF64}=ComplexF64[]
                                , method::String="nfft")
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps,:,acqData.numCoils),1)
  return [ diagOp( [ft[i] for k=1:acqData.numCoils]... )*S for i=1:acqData.numEchoes]
end

#############################
# multi echo fourier encoding
#############################
function encodingOp2d_multiEcho(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                ; slice::Int64=1
                                , correctionMap::Array{ComplexF64}=ComplexF64[]
                                , method::String="nfft")
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape, slice=slice, correctionMap=correctionMap)
  return diagOp(ft...)
end

function encodingOp3d_multiEcho(acqData::AcquisitionData, shape::NTuple{3,Int64}; correctionMap::Array{ComplexF64,3}=ComplexF64[])
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape, correctionMap=correctionMap)
  return diagOp(ft...)
end

#####################################################
# multi echo fourier encoding with coil sensitivities
#####################################################
function encodingOp2d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{2,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; slice=1
                                          , correctionMap::Array{ComplexF64}=ComplexF64[]
                                          , method::String="nfft")
  # fourier operators
  ft = encodingOps2d_simple(acqData, shape, slice=slice, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps[:,:,slice,:],:,acqData.numCoils),1)
  return diagOp( repeat(ft, acqData.numCoils)... )*S
end

function encodingOp3d_multiEcho_parallel(acqData::AcquisitionData, shape::NTuple{3,Int64}
                                          , senseMaps::Array{ComplexF64}
                                          ; correctionMap::Array{ComplexF64}=ComplexF64[]
                                          , method::String="nfft")
  # fourier operators
  ft = encodingOps3d_simple(acqData, shape, correctionMap=correctionMap)
  S = SensitivityOp(reshape(senseMaps,:,acqData.numCoils),1)
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
# function fourierEncodingOp2d(tr::Trajectory, params; subsampleIdx::Vector{Int64}=Int64[], slice::Int64=0, method::String="nfft")
#   shape = params[:shape]
#   opName="fast"
#   if get(params,:explicit,false)
#     opName="explicit"
#   end
#
#   return fourierEncodingOp2d(shape,tr,opName;subsampleIdx=subsampleIdx,slice=slice,params...)
# end

"""
return 2d Fourier encoding operator (either Explicit or NFFT)
  opname : "explicit" or "fast"
  slice : slice to which the operator will be applied
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

# function fourierEncodingOp3d(shape, tr, params; subsampleIdx::Vector{Int64}=Int64[])
#   opName="fast"
#   if get(params,:explicit,false)
#     opName="explicit"
#   end
#
#   return fourierEncodingOp3d(shape,tr,opName;params...)
# end

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
