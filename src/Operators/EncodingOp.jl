export EncodingOp, lrEncodingOp

#
# Encoding operator for one slice
#
function EncodingOp(aqData::AcquisitionData, params::Dict; set::Int64=1, numEchoes::Int64=1, parallel::Bool=false, multiEcho::Bool=false)
  # Fourier Operators
  ft = Array{AbstractLinearOperator,1}(undef,numEchoes)
  for i = 1:numEchoes
    tr = trajectory(aqData.seq,i)
    ft[i] = fourierEncodingOp(tr, params)
  end

  # coil sensitivities in case of SENSE-reconstruction
  if parallel && multiEcho
    S = SensitivityOp(params[:senseMaps][:,:,set],numEchoes)
    E = diagOp( repeat(ft, aqData.numCoils)... )*S
  elseif parallel && !multiEcho
    S = SensitivityOp(params[:senseMaps][:,:,set],1)
    E = [ diagOp( [ft[i] for k=1:aqData.numCoils]... )*S for i=1:numEchoes]
  elseif !parallel && multiEcho
    E = diagOp(ft...)
  else
    E = ft
  end

  # sampling operator in case of undersampled cartesian acquisitions
  if aqData.idx != collect(1:length(aqData.kdata))
    if get(params,:sampling, nothing) == "binary"
      M = SamplingOp(Array{Bool})(aqData.idx)
    else
      numNodes = prod(params[:shape])
      M = SamplingOp(aqData.idx,(numNodes, aqData.numEchoes, aqData.numCoils))
    end
    return M*E
  end

  return E
end


function lrEncodingOp(aqData::AcquisitionData,params::Dict; set::Int64=1, numEchoes::Int64=1, parallel::Bool=false,)

  # low rank operator
  N=prod(params[:shape])
  subspace = get( params, :phi, Matrix{Float64}(I,N,N) )
  K = size(subspace,2)
  Φ = MapSliceOp(subspace[:,:,1],2,(N, K, aqData.numCoils), (N, numEchoes, aqData.numCoils))

  # Fourier Operator
  tr = trajectory(aqData.seq,1)
  ft = fourierEncodingOp(tr, params)

  # coil sensitivities in case of SENSE-reconstruction
  if parallel
    S = SensitivityOp(params[:senseMaps][:,:,set],K)
    E = diagOp( [ft for i=1:aqData.numCoils*K]... )*S
  else
    E = diagOp( [ft for i=1:K]... )
  end

  # sampling operator in case of undersampled cartesian acquisitions
  if get(params,:sampling, nothing) == "binary"
    M = SamplingOp(Array{Bool})(aqData.idx)
  else
    M = SamplingOp(aqData.idx,(N, numEchoes, aqData.numCoils))
  end

  return M*Φ*E

end


#
# return Fourier encoding operator with(out) correction when cmap is (not) specified
#
function fourierEncodingOp(tr, params)
  densityWeighting = get(params,:densityWeighting,true)
  fft = get(params,:fft,false)
  if fft
    FFTOp(ComplexF64, params[:shape])
  elseif !haskey(params,:cmap)
    return NFFTOp(params[:shape], tr, symmetrize=densityWeighting)
  else
    echoImage = get(params, :echoImage, true)
    method = get(params, :method, "nfft")
    alpha = get(params, :alpha, 1.75)
    m = get(params, :m, 1.75)
    K = get(params, :K, 20)
    return FieldmapNFFTOp( params[:shape],  tr, params[:cmap], symmetrize=densityWeighting, echoImage=echoImage, K=K, alpha=alpha, m=m, method=method )
  end
end
