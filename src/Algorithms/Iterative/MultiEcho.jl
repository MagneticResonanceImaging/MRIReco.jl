export MultiEchoReconstruction
export AbstractMultiEchoParameters
export MultiEchoIterativeParameters

abstract type AbstractMultiEchoParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct MultiEchoReconstruction{P<:AbstractMultiEchoParameters, C <: AbstractIterativeMRIRecoContextParameter{P}} <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::C
end

function (weighting::AbstractMRIRecoWeightingParameters)(algo::MultiEchoReconstruction)
  weights = weighting(AbstractIterativeMRIRecoAlgorithm)
  return vcat(weights...)
end


function (sparsity::AbstractSparsityParameters)(algo::MultiEchoReconstruction, numTerms::Int64)
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

function (ep::EncodingParameters)(::MultiEchoReconstruction, slice::Int)
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  
  encParams = buildEncodingParams(ep, S)
  return encodingOp_multiEcho(acqData, reconSize; slice=slice, encParams...)
end


"""
    MultiEchoIterativeParameters{E, W, S} <: AbstractMultiEchoParameters

Parameters for multi-echo iterative MRI reconstruction.
Reconstructs all echoes jointly using a combined operator.

# Type Parameters
- `E` - Encoding parameters type
- `W` - Weighting parameters type  
- `S` - Solver parameters type

# Fields
- `encodingParams::E` - Encoding operator parameters
- `weightingParams::W` - Sampling weighting parameters
- `solverParams::S` - Least squares solver parameters

# Callable Interface

    # Allocation - allocates output array and returns iteration indices and weights
    (params::MultiEchoIterativeParameters)(algo, reconSize) -> (Ireco, indices, weights)

    # Loop body - called per index with weights and output array
    (params::MultiEchoIterativeParameters)(algo, Ireco, index, weights)

    # Finalization - reshapes and wraps result
    (params::MultiEchoIterativeParameters)(algo, Ireco) -> AxisArray

# Usage
```julia
# Create the algorithm parameter
algoParams = MultiEchoIterativeParameters(
  encoding = EncodingParameters(),
  weighting = DensityWeightingParameters(),
  solver = LeastSquaresSolverParameter(solver=FISTA)
)

# Wrap in context (serial or threaded)
reco = SerialIterativeMRIRecoContextParameter(; parameter=algoParams)
# or
reco = ThreadedIterativeMRIRecoContextParameter(; parameter=algoParams, scheduler=DynamicScheduler())
```
"""
@parameter struct MultiEchoIterativeParameters{
    E <: AbstractMRIRecoEncodingParameters,
    W <: AbstractMRIRecoWeightingParameters,
    S <: LeastSquaresSolverParameter
} <: AbstractMultiEchoParameters
  encodingParams::E
  weightingParams::W
  solverParams::S
end

function (params::MultiEchoIterativeParameters)(
  ::MultiEchoReconstruction, 
  reconSize::NTuple{D, Int64}
) where D
  acqData = ctx_acqData()
  T = real(eltype(ctx_storageType()))
  
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  numSl = numSlices(acqData)
  numRep = numRepetitions(acqData)
  
  Ireco = zeros(Complex{T}, prod(reconSize) * numContr, numChan, numSl, numRep)
  indices = CartesianIndices((numRep, numSl))
  
  weights = params.weightingParams(MultiEchoReconstruction)
  
  return (Ireco, indices, weights)
end

function (params::MultiEchoIterativeParameters)(
  algoT::MultiEchoReconstruction, 
  Ireco::Array{Complex{T}, 4},
  index::CartesianIndex{2}, 
  weights::AbstractVector
) where T
  reconSize = ctx_reconSize()
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()
  
  l, i = index[1], index[2]  # rep, slice
  
  numChan = numChannels(acqData)
  
  F = params.encodingParams(algoT, i)
  
  W = WeightingOp(Complex{T}; weights=weights)
  
  for j in 1:numChan
    kdata = arrayType(multiEchoData(acqData, j, i, rep=l)) .* weights
    EFull = ProdOp(W, F)
    EFullᴴEFull = normalOperator(EFull; normalOpParams(arrayType)...)
    
    I = params.solverParams(algoT, kdata, EFull, EFullᴴEFull)
    
    Ireco[:, j, i, l] = Array(I)
  end
end

function (params::MultiEchoIterativeParameters)(
  algo::MultiEchoReconstruction, 
  Ireco::Array{Complex{T}, 4}
) where T
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  
  numSl = numSlices(acqData)
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  numRep = numRepetitions(acqData)
  
  encDims = ndims(trajectory(acqData))
  
  if encDims == 2
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], numContr, numChan, numSl, numRep)
    Ireco_ = permutedims(Ireco_, [1, 2, 5, 3, 4, 6])
  else
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numContr, numChan, numRep)
  end
  
  return makeAxisArray(Ireco_, acqData)
end