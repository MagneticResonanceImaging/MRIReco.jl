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
- `algo` - The algorithm
- `reconSize::NTuple{D, Int64}` - Reconstruction size
- Returns a tuple of `(Ireco, indices, weights, extra...)` where:
  - `Ireco::Array{Complex{T}}` - Output image array
  - `indices` - CartesianIndices for iteration
  - `weights` - Sampling weights
  - `extra...` - Optional additional data (e.g., decorrelatedSenseMaps, L_inv for coil recon)

## Loop Body
    (params)(algo, Ireco, index, weights, extra...)

Called once per iteration index. Writes results to Ireco.
- `algo` - The algorithm
- `Ireco` - Output array to write to
- `index::CartesianIndex` - Current index, such as (rep, slice)
- `weights` - Sampling weights from allocation
- `extra...` - Optional additional data from allocation

## Finalization
    (params)(algo, Ireco) -> result

Finalizes and returns the reconstruction result.
- `algo` - The algorithm
- `Ireco::Array{Complex{T}}` - The filled output array
- Returns final result (typically AxisArray)
"""
abstract type AbstractIterativeRecoParameters <: AbstractMRIRecoParameters end

include("IterativeContextParameters.jl")

include("Standard.jl")
include("MultiEcho.jl")
include("MultiCoil.jl")
include("MultiCoilMultiEcho.jl")
include("MultiCoilMultiEchoSubspace.jl")
