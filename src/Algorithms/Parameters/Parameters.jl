include("SolverParameters.jl")
include("WeightingParameters.jl")
include("EncodingParameters.jl")

export AbstractIterativeRecoParameters

"""
    AbstractIterativeRecoParameters <: AbstractMRIRecoParameters

Abstract base type for iterative reconstruction algorithm parameters.
Subtypes must implement callable methods for use with 
`AbstractIterativeMRIRecoContextParameter`.

# Required Callable Interface

## Allocation
    (params)(algo, reconSize) -> (Ireco, indices)

Allocates the output array and returns iteration indices.
- `algo` - Algorithm (type)
- `reconSize::NTuple{D, Int64}` - Reconstruction size
- Returns a tuple of `(Ireco, indices)` where:
  - `Ireco::Array{Complex{T}, 5}` - Output image array
  - `indices` - CartesianIndices for iteration

## Loop Body
    (params)(algo, Ireco, index, weights)
    (params)(algo, Ireco, index)

Called once per iteration index. Writes results to Ireco.
- `algo` - Algorithm (type)
- `Ireco` - Output array to write to
- `index::CartesianIndex` - Current index (rep, slice)
- `weights` - Pre-computed weights from weighting parameter, optional if getWeighting is provided

## Finalization
    (params)(algo, Ireco) -> result

Finalizes and returns the reconstruction result.
- `algo` - Algorithm (type)
- `Ireco::Array{Complex{T}, 5}` - The filled output array
- Returns final result (typically AxisArray)

# Optional Accessor

    getWeighting(params) -> weighting::AbstractMRIRecoWeightingParameters

Returns the weighting parameter used to compute weights for the loop-body.
The context parameter calls `getWeighting(params)` to get actual weights.
"""
abstract type AbstractIterativeRecoParameters <: AbstractMRIRecoParameters end

getWeighting(params::AbstractIterativeRecoParameters) = nothing

include("IterativeContextParameters.jl")

