export EncodingOp, lrEncodingOp

#
# Encoding operator for one slice
#
function EncodingOp(acqData::AcquisitionData, params::Dict; slice=1,
                           parallel=false, multiEcho=false)
  # Fourier Operators
  ft = Array{AbstractLinearOperator,1}(undef,acqData.numEchoes)

  for i = 1:acqData.numEchoes
    tr = trajectory(acqData.seq,i)
    ft[i] = fourierEncodingOp(tr, params, slice=slice)
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
  if acqData.idx != collect(1:length(acqData.kdata))
    if get(params,:sampling, nothing) == "binary"
      M = SamplingOp(Array{Bool})(acqData.idx)
    else
      numNodes = prod(params[:shape])
      M = SamplingOp(acqData.idx,(numNodes, acqData.numEchoes, acqData.numCoils))
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
  tr = trajectory(acqData.seq,1)
  ft = fourierEncodingOp(tr, params)

  # coil sensitivities in case of SENSE-reconstruction
  if parallel
    S = SensitivityOp(params[:senseMaps][:,:],K)
    E = diagOp( [ft for i=1:acqData.numCoils*K]... )*S
  else
    E = diagOp( [ft for i=1:K]... )
  end

  # sampling operator in case of undersampled cartesian acquisitions
  if get(params,:sampling, nothing) == "binary"
    M = SamplingOp(Array{Bool})(acqData.idx)
  else
    M = SamplingOp(acqData.idx,(N, numEchoes, acqData.numCoils))
  end

  return M*Φ*E

end


#
# return Fourier encoding operator with(out) correction when cmap is (not) specified
#
function fourierEncodingOp(tr, params; slice=1)
  densityWeighting = get(params,:densityWeighting,true)
  fft = get(params,:fft,false)
  if fft
    FFTOp(ComplexF64, params[:shape])
  elseif !haskey(params,:correctionMap)
    return NFFTOp(params[:shape], tr, symmetrize=densityWeighting)
  else
    echoImage = get(params, :echoImage, true)
    method = get(params, :method, "nfft")
    alpha = get(params, :alpha, 1.75)
    m = get(params, :m, 1.75)
    K = get(params, :K, 20)
    return FieldmapNFFTOp( params[:shape],  tr, params[:correctionMap][:,:,slice],
                           symmetrize=densityWeighting, echoImage=echoImage, K=K,
                           alpha=alpha, m=m, method=method )
  end
end
