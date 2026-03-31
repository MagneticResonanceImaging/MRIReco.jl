export MultiCoilReconstruction

export AbstractMultiCoilParameters
export MultiCoilIterativeParameters

abstract type AbstractMultiCoilParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct MultiCoilReconstruction{P<:AbstractMultiCoilParameters, C <: AbstractIterativeMRIRecoContextParameter{P}} <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::C
end

function (ep::EncodingParameters)(::Type{<:MultiCoilReconstruction}, slice::Int, decorrelatedSenseMaps)
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  
  encParams = buildEncodingParams(ep, S)
  return encodingOps_parallel(acqData, reconSize, decorrelatedSenseMaps; slice=slice, encParams...)
end


"""
    MultiCoilIterativeParameters{E, W, S, C} <: AbstractMultiCoilParameters

Parameters for multi-coil (SENSE) iterative MRI reconstruction.
Reconstructs coil images jointly using sensitivity encoding.

# Type Parameters
- `E` - Encoding parameters type
- `W` - Weighting parameters type  
- `S` - Solver parameters type
- `C` - Coil parameters type

# Fields
- `encodingParams::E` - Encoding operator parameters
- `weightingParams::W` - Sampling weighting parameters
- `solverParams::S` - Least squares solver parameters
- `coilParams::C` - Coil sensitivity and noise decorrelation parameters

# Callable Interface

    # Allocation - allocates output array and returns iteration indices, weights, and coil data
    (params::MultiCoilIterativeParameters)(algo, reconSize) -> (Ireco, indices, weights, decorrelatedSenseMaps, L_inv)

    # Loop body - called per index with weights and output array
    (params::MultiCoilIterativeParameters)(algo, Ireco, index, weights, decorrelatedSenseMaps, L_inv)

    # Finalization - reshapes and wraps result
    (params::MultiCoilIterativeParameters)(algo, Ireco) -> AxisArray

# Usage
```julia
# Create the algorithm parameter
algoParams = MultiCoilIterativeParameters(
  encoding = EncodingParameters(),
  weighting = DensityWeightingParameters(),
  solver = LeastSquaresSolverParameter(solver=FISTA),
  coilParams = CoilParameters(; senseMaps=smaps, noiseData=noiseData)
)

# Wrap in context (serial or threaded)
reco = SerialIterativeMRIRecoContextParameter(; parameter=algoParams)
# or
reco = ThreadedIterativeMRIRecoContextParameter(; parameter=algoParams, scheduler=DynamicScheduler())
```
"""
@parameter struct MultiCoilIterativeParameters{
    E <: AbstractMRIRecoEncodingParameters,
    W <: AbstractMRIRecoWeightingParameters,
    S <: LeastSquaresSolverParameter,
    C <: AbstractCoilParameters
} <: AbstractMultiCoilParameters
  encodingParams::E
  weightingParams::W
  solverParams::S
  coilParams::C
end

function (params::MultiCoilIterativeParameters)(
  algo::MultiCoilReconstruction, 
  reconSize::NTuple{D, Int64}
) where D
  acqData = ctx_acqData()
  T = real(eltype(ctx_storageType()))
  
  numContr = numContrasts(acqData)
  numSl = numSlices(acqData)
  numRep = numRepetitions(acqData)
  
  Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numRep)
  indices = CartesianIndices((numRep, numSl, numContr))
  
  L_inv, decorrelatedSenseMaps = params.coilParams(algo)
  weights = params.weightingParams(algo)
  
  return (Ireco, indices, weights, decorrelatedSenseMaps, L_inv)
end

function (params::MultiCoilIterativeParameters)(
  algoT::Type{<:MultiCoilReconstruction},
  Ireco::Array{Complex{T}, 4},
  index::CartesianIndex{3}, 
  weights::AbstractVector,
  decorrelatedSenseMaps,
  L_inv
) where T
  reconSize = ctx_reconSize()
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()
  
  l, k, j = index[1], index[2], index[3]  # rep, slice, contrast
  
  numChan = numChannels(acqData)
  
  E = params.encodingParams(algoT, k, decorrelatedSenseMaps)
  
  W = WeightingOp(Complex{T}; weights=weights[j], rep=numChan)
  kdata = arrayType(multiCoilData(acqData, j, k, rep=l)) .* repeat(weights[j], numChan)
  
  if !isnothing(L_inv)
    kdata = vec(reshape(kdata, :, numChan) * L_inv')
  end
  
  EFull = ∘(W, E[j])
  EFullᴴEFull = normalOperator(EFull; normalOpParams(arrayType)...)
  
  I = params.solverParams(algoT, kdata, EFull, EFullᴴEFull)
  
  if isCircular(trajectory(acqData, j))
    circularShutter!(reshape(I, reconSize), 1.0)
  end
  
  Ireco[:, k, j, l] = Array(I)
end

function (params::MultiCoilIterativeParameters)(
  algo::MultiCoilReconstruction, 
  Ireco::Array{Complex{T}, 4}
) where T
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  
  numSl = numSlices(acqData)
  numContr = numContrasts(acqData)
  numRep = numRepetitions(acqData)
  
  Ireco_ = reshape(Ireco, volumeSize(reconSize, numSl)..., numContr, 1, numRep)
  return makeAxisArray(Ireco_, acqData)
end