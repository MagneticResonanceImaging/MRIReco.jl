export AbstractIterativeMRIRecoContextParameter
export SerialIterativeMRIRecoContextParameter, ThreadedIterativeMRIRecoContextParameter

"""
    AbstractIterativeMRIRecoContextParameter{P} <: AbstractMRIRecoParameters

Abstract base type for iterative reconstruction context parameters.
Provides the runtime framework (context, parallelization, allocation) 
while the inner parameter implements the algorithm-specific loop body.

# Type Parameters
- `P` - The inner algorithm parameter type (e.g., AbstractStandardParameters)
"""
abstract type AbstractIterativeMRIRecoContextParameter{P <: AbstractIterativeRecoParameters} <: AbstractMRIRecoParameters end

"""
    SerialIterativeMRIRecoContextParameter{P, RS, AT} <: AbstractIterativeMRIRecoContextParameter{P}

Serial (non-parallel) iterative reconstruction context parameter.
Uses a simple for loop for iteration.

# Type Parameters
- `P` - The inner algorithm parameter type
- `RS` - ReconSize type (NTuple or Nothing)
- `AT` - Array type

# Fields
- `reconSize::RS` - Reconstruction size (nothing = infer from acqData)
- `arrayType::Type{AT}` - Array type for computation
- `parameter::P` - Inner algorithm-specific parameter
"""
@parameter struct SerialIterativeMRIRecoContextParameter{
    P <: AbstractIterativeRecoParameters,
    RS <: Union{Nothing, NTuple{D, Int64} where D},
    AT <: AbstractArray
} <: AbstractIterativeMRIRecoContextParameter{P}
  reconSize::RS = nothing
  arrayType::Type{AT} = Array
  parameter::P
end

"""
    ThreadedIterativeMRIRecoContextParameter{P, RS, AT, S} <: AbstractIterativeMRIRecoContextParameter{P}

Threaded iterative reconstruction context parameter using OhMyThreads.
Uses `@tasks` for parallel execution.

# Type Parameters
- `P` - The inner algorithm parameter type
- `RS` - ReconSize type (NTuple or Nothing)
- `AT` - Array type
- `S` - Scheduler type

# Fields
- `reconSize::RS` - Reconstruction size (nothing = infer from acqData)
- `arrayType::Type{AT}` - Array type for computation
- `scheduler::S` - Task scheduler (nothing = auto-based on arrayType)
- `parameter::P` - Inner algorithm-specific parameter
"""
@parameter struct ThreadedIterativeMRIRecoContextParameter{
    P <: AbstractIterativeRecoParameters,
    RS <: Union{Nothing, NTuple{D, Int64} where D},
    AT <: AbstractArray,
    S <: Union{Nothing, Symbol, OhMyThreads.Scheduler}
} <: AbstractIterativeMRIRecoContextParameter{P}
  reconSize::RS = nothing
  arrayType::Type{AT} = Array
  scheduler::S = nothing
  parameter::P
end

function setupIterativeReconSize(acqData::AcquisitionData, reconSize::Nothing)
  rs = encodingSize(acqData)
  return setupIterativeReconSize(acqData, rs)
end

function setupIterativeReconSize(acqData::AcquisitionData, reconSize::NTuple)
  encDims = ndims(trajectory(acqData,1))
  
  if encDims == 3 && numSlices(acqData) > 1
    error("reconstruction of multiple 3d-encoded volumina is not yet supported")
  end
  
  red3d = encDims == 2 && length(reconSize) == 3
  if red3d
    reconSize = (reconSize[2], reconSize[3])
  end
  
  return reconSize
end

get_scheduler(::Nothing, arrayType) = scheduler(arrayType)
get_scheduler(scheduler, arrayType) = scheduler

# Serial implementation
function (params::SerialIterativeMRIRecoContextParameter{P})(algo::AbstractIterativeMRIRecoAlgorithm, acqData::AcquisitionData) where P <: AbstractIterativeRecoParameters
  reconSize = setupIterativeReconSize(acqData, params.reconSize)
  
  encDims = ndims(trajectory(acqData))
  if encDims != length(reconSize)
    error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
  end
  
  with(MRIRECO_CONTEXT => MRIRecoContext(reconSize, acqData, params.arrayType)) do
    # Allocation - call algorithm parameter (returns Ireco, indices, weights, extra...)
    Ireco, indices, extra... = params.parameter(algo, reconSize)
    
    # Loop over indices
    for index in indices
      params.parameter(algo, Ireco, index, extra...)
    end
    
    # Finalization - call algorithm parameter
    return params.parameter(algo, Ireco)
  end
end

# Threaded implementation
function (params::ThreadedIterativeMRIRecoContextParameter{P})(algo::AbstractIterativeMRIRecoAlgorithm, acqData::AcquisitionData) where P <: AbstractIterativeRecoParameters
  reconSize = setupIterativeReconSize(acqData, params.reconSize)
  
  encDims = ndims(trajectory(acqData))
  if encDims != length(reconSize)
    error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
  end
  
  sched = get_scheduler(params.scheduler, params.arrayType)
  
  with(MRIRECO_CONTEXT => MRIRecoContext(reconSize, acqData, params.arrayType)) do
    # Allocation - call algorithm parameter (returns Ireco, indices, weights, extra...)
    Ireco, indices, weights, extra... = params.parameter(algo, reconSize)
    
    # Loop over indices with @tasks
    @tasks for index in indices
      @set scheduler = sched
      params.parameter(algo, Ireco, index, weights, extra...)
    end
    
    # Finalization - call algorithm parameter
    return params.parameter(algo, Ireco)
  end
end