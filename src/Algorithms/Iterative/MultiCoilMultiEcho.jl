export MultiCoilMultiEchoReconstruction

@reconstruction mutable struct MultiCoilMultiEchoReconstruction{P<:AbstractMultiCoilParameters, C <: AbstractIterativeMRIRecoContextParameter{P}} <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::C
end

function (weighting::AbstractMRIRecoWeightingParameters)(algo::MultiCoilMultiEchoReconstruction)
  numChan = numChannels(ctx_acqData())
  weights = weighting(AbstractIterativeMRIRecoAlgorithm)
  return repeat(vcat(weights...), numChan)
end

function (sparsity::AbstractSparsityParameters)(algo::MultiCoilMultiEchoReconstruction, numTerms::Int64)
  trafos = sparsity(AbstractIterativeMRIRecoAlgorithm, numTerms)
  acqData = ctx_acqData()
  numContr = numContrasts(acqData)
  
  result = []
  for trafo in trafos
    if isnothing(trafo)
      push!(result, nothing)
    else
      push!(result, DiagOp(repeat([trafo], numContr)...))
    end
  end
  return result
end

function (ep::EncodingParameters)(::Type{<:MultiCoilMultiEchoReconstruction}, slice::Int, decorrelatedSenseMaps)
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  
  encParams = buildEncodingParams(ep, S)
  return encodingOp_multiEcho_parallel(acqData, reconSize, decorrelatedSenseMaps; slice=slice, encParams...)
end

function (params::MultiCoilIterativeParameters)(
  algo::MultiCoilMultiEchoReconstruction, 
  reconSize::NTuple{D, Int64}
) where D
  acqData = ctx_acqData()
  T = real(eltype(ctx_storageType()))
  
  numContr = numContrasts(acqData)
  numSl = numSlices(acqData)
  numRep = numRepetitions(acqData)
  
  Ireco = zeros(Complex{T}, prod(reconSize) * numContr, numSl, numRep)
  indices = CartesianIndices((numRep, numSl))
  
  L_inv, decorrelatedSenseMaps = params.coilParams(algo)
  weights = params.weightingParams(algo)
  
  return (Ireco, indices, weights, decorrelatedSenseMaps, L_inv)
end

function (params::MultiCoilIterativeParameters)(
  algoT::Type{<:MultiCoilMultiEchoReconstruction},
  Ireco::Array{Complex{T}, 3},
  index::CartesianIndex{2}, 
  weights::AbstractVector{Complex{T}},
  decorrelatedSenseMaps,
  L_inv
) where T
  reconSize = ctx_reconSize()
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()
  
  l, i = index[1], index[2]  # rep, slice
  
  numChan = numChannels(acqData)
  
  E = params.encodingParams(algoT, i, decorrelatedSenseMaps)
  
  W = WeightingOp(Complex{T}; weights=weights)
  kdata = arrayType(multiCoilMultiEchoData(acqData, i)) .* weights
  
  if !isnothing(L_inv)
    kdata = vec(reshape(kdata, :, numChan) * L_inv')
  end
  
  EFull = ∘(W, E)
  EFullᴴEFull = normalOperator(EFull; normalOpParams(arrayType)...)
  
  I = params.solverParams(algoT, kdata, EFull, EFullᴴEFull)
  
  Ireco[:, i, l] = Array(I)
end

function (params::MultiCoilIterativeParameters)(
  algo::MultiCoilMultiEchoReconstruction, 
  Ireco::Array{Complex{T}, 3}
) where T
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  
  numSl = numSlices(acqData)
  numContr = numContrasts(acqData)
  numRep = numRepetitions(acqData)
  
  encDims = ndims(trajectory(acqData))
  
  if encDims == 2
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], numContr, numSl, 1, numRep)
    Ireco_ = permutedims(Ireco_, [1, 2, 4, 3, 5, 6])
  else
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numContr, 1, numRep)
  end
  
  return makeAxisArray(Ireco_, acqData)
end