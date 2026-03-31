include("SolverParameters.jl")
include("WeightingParameters.jl")
include("EncodingParameters.jl")
include("CoilParameters.jl")

export AbstractIterativeRecoParameters

"""
    AbstractIterativeRecoParameters <: AbstractMRIRecoParameters

Abstract base type for iterative reconstruction algorithm parameters.
Subtypes must implement callable methods for use with 
`AbstractIterativeMRIRecoContextParameter`.

# Required Callable Interface

## Allocation
    (params)(algo, reconSize) -> (Ireco, indices, weights, extra...)

Allocates the output array and returns iteration indices and weights.
- `algo` - Algorithm (type)
- `reconSize::NTuple{D, Int64}` - Reconstruction size
- Returns a tuple of `(Ireco, indices, extra...)` where:
  - `Ireco::Array{Complex{T}}` - Output image array
  - `indices` - CartesianIndices for iteration
  - `extra...` - Optional additional data (e.g., decorrelatedSenseMaps, L_inv for coil recon)

## Loop Body
    (params)(algo, Ireco, index, extra...)

Called once per iteration index. Writes results to Ireco.
- `algo` - Algorithm (type)
- `Ireco` - Output array to write to
- `index::CartesianIndex` - Current index (rep, slice)
- `extra...` - Optional additional data from allocation

## Finalization
    (params)(algo, Ireco) -> result

Finalizes and returns the reconstruction result.
- `algo` - Algorithm (type)
- `Ireco::Array{Complex{T}}` - The filled output array
- Returns final result (typically AxisArray)
"""
abstract type AbstractIterativeRecoParameters <: AbstractMRIRecoParameters end

include("IterativeContextParameters.jl")

