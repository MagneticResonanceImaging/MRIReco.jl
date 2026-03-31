export StandardReconstruction
export StandardIterativeParameters

export AbstractStandardParameters

abstract type AbstractStandardParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct StandardReconstruction{P<:AbstractStandardParameters, C <: AbstractIterativeMRIRecoContextParameter{P}} <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::C
end

"""
    StandardIterativeParameters{E, W, S} <: AbstractStandardParameters

Parameters for standard iterative MRI reconstruction.
Only implements the algorithm-specific loop body - context is handled by 
SerialIterativeMRIRecoContextParameter or ThreadedIterativeMRIRecoContextParameter.

# Type Parameters
- `E` - Encoding parameters type
- `W` - Weighting parameters type  
- `S` - Solver parameters type

# Fields
- `encoding::E` - Encoding operator parameters
- `weighting::W` - Sampling weighting parameters
- `solver::S` - Least squares solver parameters

# Callable Interface

    # Allocation - allocates output array and returns iteration indices and weights
    (params::StandardIterativeParameters)(algo, reconSize) -> (Ireco, indices, weights)

    # Loop body - called per index with weights and output array
    (params::StandardIterativeParameters)(algo, Ireco, index, weights)

    # Finalization - reshapes and wraps result
    (params::StandardIterativeParameters)(algo, Ireco) -> AxisArray

# Usage
```julia
# Create the algorithm parameter
algoParams = StandardIterativeParameters(
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
@parameter struct StandardIterativeParameters{
    E <: AbstractMRIRecoEncodingParameters,
    W <: AbstractMRIRecoWeightingParameters,
    S <: LeastSquaresSolverParameter
} <: AbstractStandardParameters
  encodingParams::E
  weightingParams::W
  solverParams::S
end

function (params::StandardIterativeParameters)(
  algoT::Type{<:StandardReconstruction}, 
  reconSize::NTuple{D, Int64}
) where D
  acqData = ctx_acqData()
  T = real(eltype(ctx_storageType()))
  
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  numSl = numSlices(acqData)
  numRep = numRepetitions(acqData)
  
  Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numChan, numRep)
  indices = CartesianIndices((numRep, numSl))
  
  weights = params.weightingParams(algoT)
  
  return (Ireco, indices, weights)
end

function (params::StandardIterativeParameters)(
  algoT::Type{<:StandardReconstruction}, 
  Ireco::Array{Complex{T}, 5},
  index::CartesianIndex{2}, 
  weights
) where T
  reconSize = ctx_reconSize()
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()
  
  l, k = index[1], index[2]  # rep, slice
  
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  
  F = params.encodingParams(algoT, k)
  
  for j in 1:numContr
    W = WeightingOp(Complex{T}; weights=weights[j])
    for i in 1:numChan
      kdata = arrayType(kData(acqData, j, i, k, rep=l)) .* weights[j]
      EFull = ProdOp(W, F[j])
      EFullᴴEFull = normalOperator(EFull; normalOpParams(arrayType)...)
      
      I = params.solverParams(algoT, kdata, EFull, EFullᴴEFull)
      
      if isCircular(trajectory(acqData, j))
        circularShutter!(reshape(I, reconSize), 1.0)
      end
      
      Ireco[:, k, j, i, l] = Array(I)
    end
  end
end

function (params::StandardIterativeParameters)(
  algoT::Type{<:StandardReconstruction}, 
  Ireco::Array{Complex{T}, 5}
) where T
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  
  numSl = numSlices(acqData)
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  numRep = numRepetitions(acqData)
  
  Ireco = reshape(Ireco, volumeSize(reconSize, numSl)..., numContr, numChan, numRep)
  return makeAxisArray(Ireco, acqData)
end