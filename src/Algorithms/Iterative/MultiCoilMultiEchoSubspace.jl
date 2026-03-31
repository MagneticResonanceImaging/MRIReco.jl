export MultiCoilMultiEchoSubspaceReconstruction

export AbstractSubspaceParameters
abstract type AbstractSubspaceParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct MultiCoilMultiEchoSubspaceReconstruction{
    P<:AbstractMultiCoilParameters, 
    C<:AbstractIterativeMRIRecoContextParameter{P}
} <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::C
  numBasis::Union{Nothing, Int64} = nothing

  @init function getBasisNumber(algo::MultiCoilMultiEchoSubspaceReconstruction)
    encoding = algo.parameter.parameter.encodingParams
    if encoding isa SubspaceEncodingParameters
      algo.numBasis = size(encoding.basis, 2)
    else
      throw(ArgumentError("Subspace reconstruction requires SubspaceEncodingParameters, found $(typeof(encoding))"))
    end
  end
end

function (weighting::AbstractMRIRecoWeightingParameters)(algo::MultiCoilMultiEchoSubspaceReconstruction)
  weights = weighting(AbstractIterativeMRIRecoAlgorithm)
  return vcat(weights...)
end

function (sparsity::AbstractSparsityParameters)(algo::MultiCoilMultiEchoSubspaceReconstruction, numTerms::Int64)
  trafos = sparsity(AbstractIterativeMRIRecoAlgorithm, numTerms)
  acqData = ctx_acqData()
  numContr = numContrasts(acqData)
  
  numBasis = algo.numBasis
  
  result = []
  for trafo in trafos
    if isnothing(trafo)
      push!(result, nothing)
    else
      push!(result, DiagOp(repeat([trafo], numBasis)...))
    end
  end
  return result
end

function (ep::EncodingParameters)(::Type{<:MultiCoilMultiEchoSubspaceReconstruction}, slice::Int, decorrelatedSenseMaps)
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  
  encParams = buildEncodingParams(ep, S)
  return encodingOp_multiEcho_parallel(acqData, reconSize, decorrelatedSenseMaps; slice=slice, encParams...)
end

function (params::MultiCoilIterativeParameters)(
  algo::MultiCoilMultiEchoSubspaceReconstruction, 
  reconSize::NTuple{D, Int64}
) where D
  acqData = ctx_acqData()
  T = real(eltype(ctx_storageType()))
  
  numContr = numContrasts(acqData)
  numSl = numSlices(acqData)
  numRep = numRepetitions(acqData)
  
  numBasis = algo.numBasis
  
  Ireco = zeros(Complex{T}, prod(reconSize) * numBasis, numSl, numRep)
  indices = CartesianIndices((numRep, numSl))
  
  L_inv, decorrelatedSenseMaps = params.coilParams(algo)
  weights = params.weightingParams(algo)
  
  return (Ireco, indices, weights, decorrelatedSenseMaps, L_inv)
end

function (params::MultiCoilIterativeParameters)(
  algoT::Type{<:MultiCoilMultiEchoSubspaceReconstruction},
  Ireco::Array{Complex{T}, 3},
  index::CartesianIndex{2}, 
  weights::AbstractVector,
  decorrelatedSenseMaps,
  L_inv
) where T
  reconSize = ctx_reconSize()
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()
  
  l, i = index[1], index[2]  # rep, slice
  
  numChan = numChannels(acqData)
  
  E = params.encodingParams(algoT, i, decorrelatedSenseMaps)
  
  W = WeightingOp(Complex{T}; weights=weights, rep=numChan)
  kdata = arrayType(multiCoilMultiEchoData(acqData, i)) .* repeat(weights, numChan)
  
  if !isnothing(L_inv)
    kdata = vec(reshape(kdata, :, numChan) * L_inv')
  end
  
  EFull = ∘(W, E)
  EFullᴴEFull = normalOperator(EFull; normalOpParams(arrayType)...)
  
  I = params.solverParams(algoT, kdata, EFull, EFullᴴEFull)
  
  Ireco[:, i, l] = Array(I)
end

function (params::MultiCoilIterativeParameters)(
  algo::MultiCoilMultiEchoSubspaceReconstruction, 
  Ireco::Array{Complex{T}, 3}
) where T
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  
  numSl = numSlices(acqData)
  numContr = numContrasts(acqData)
  numRep = numRepetitions(acqData)
  
  numBasis = algo.numBasis
  
  encDims = ndims(trajectory(acqData))
  
  if encDims == 2
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], numBasis, numSl, 1, numRep)
    Ireco_ = permutedims(Ireco_, [1, 2, 4, 3, 5, 6])
  else
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numBasis, 1, numRep)
  end
  
  return makeAxisArray(Ireco_, acqData)
end