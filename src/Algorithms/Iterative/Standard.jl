export StandardReconstruction
export StandardIterativeParameters
export getWeighting

export AbstractStandardParameters

abstract type AbstractStandardParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct StandardReconstruction{P<:AbstractStandardParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
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

    # Allocation - allocates output array and returns iteration indices
    (params::StandardIterativeParameters)(algo, reconSize) -> (Ireco, indices)

    # Loop body - called per index with weights and output array
    (params::StandardIterativeParameters)(algo, index, weights, Ireco)

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

# Run reconstruction
image = reco(acqData)
```
"""
@parameter struct StandardIterativeParameters{
    E <: AbstractMRIRecoEncodingParameters,
    W <: AbstractMRIRecoWeightingParameters,
    S <: LeastSquaresSolverParameter
} <: AbstractStandardParameters
  encoding::E
  weighting::W
  solver::S
end

function StandardIterativeParameters(;
  encoding = EncodingParameters(),
  weighting = DensityWeightingParameters(),
  solver = LeastSquaresSolverParameter()
)
  StandardIterativeParameters(encoding, weighting, solver)
end

# Accessor for weighting parameter
getWeighting(params::StandardIterativeParameters) = params.weighting

# Level 1: Allocation - allocates the output array and returns indices for iteration
function (params::StandardIterativeParameters)(
  algoT::Type{<:AbstractIterativeMRIRecoAlgorithm}, 
  reconSize::NTuple{D, Int64}
) where D
  acqData = ctx_acqData()
  T = real(eltype(acqData))
  
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  numSl = numSlices(acqData)
  numRep = numRepetitions(acqData)
  
  Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numChan, numRep)
  indices = CartesianIndices((numRep, numSl))
  
  return (Ireco, indices)
end

# Level 2: Loop body - called once per index (slice/rep combination)
# The index is a CartesianIndex (rep, slice)
# Weights are pre-computed and passed in
# Ireco is the output array to write to
function (params::StandardIterativeParameters)(
  algoT::Type{<:AbstractIterativeMRIRecoAlgorithm}, 
  index::CartesianIndex{2}, 
  weights,
  Ireco::Array{Complex{T}, 5}
) where T
  reconSize = ctx_reconSize()
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()
  
  l, k = index[1], index[2]  # rep, slice
  
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  
  F = params.encoding(algoT, k)
  
  for j in 1:numContr
    W = WeightingOp(Complex{T}; weights=weights[j])
    for i in 1:numChan
      kdata = arrayType(kData(acqData, j, i, k, rep=l)) .* weights[j]
      EFull = ProdOp(W, F[j])
      EFullᴴEFull = normalOperator(EFull; normalOpParams(arrayType)...)
      
      I = params.solver(algoT, kdata, EFull, EFullᴴEFull)
      
      if isCircular(trajectory(acqData, j))
        circularShutter!(reshape(I, reconSize), 1.0)
      end
      
      Ireco[:, k, j, i, l] = Array(I)
    end
  end
end

# Level 3: Finalization - reshapes and wraps result in AxisArray
function (params::StandardIterativeParameters)(
  algoT::Type{<:AbstractIterativeMRIRecoAlgorithm}, 
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