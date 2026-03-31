export AbstractMRIRecoEncodingParameters
export EncodingParameters, CustomEncodingParameters, SubspaceEncodingParameters

"""
    AbstractMRIRecoEncodingParameters <: AbstractMRIRecoParameters

Abstract base type for encoding operator parameters used in MRI reconstruction.
"""
abstract type AbstractMRIRecoEncodingParameters <: AbstractMRIRecoParameters end

"""
    EncodingParameters{T} <: AbstractMRIRecoEncodingParameters

Parameters for computing encoding operators from acquisition data.

# Type Parameters
- `C` - The type of the correctionMap (Union{AbstractArray{Complex}, Nothing})

# Fields
- `correctionMap::C` - Fieldmap for off-resonance correction
- `method::Union{String, Nothing}` - Encoding method (typically "fast" for NFFT)
- `toeplitz::Union{Bool, Nothing}` - Enable Toeplitz
- `oversamplingFactor::Union{Float64, Nothing}` - NFFT oversampling factor
- `kernelSize::Union{Int, Nothing}` - NFFT kernel size
- `K::Union{Int, Nothing}` - Number of translates for time-segmentation
- `K_tol::Union{Float64, Nothing}` - Tolerance for LS approach

# Callable Interface
    (ep::EncodingParameters)(::Type{<:AbstractMRIRecoAlgorithm}, slice::Int)

# Returns
- Vector of encoding operators for the given slice

# Example
```julia
ep = EncodingParameters()
encOps = ep(AbstractIterativeMRIRecoAlgorithm, 1)

# With fieldmap correction
ep = EncodingParameters(; correctionMap=fieldmap)
encOps = ep(AbstractIterativeMRIRecoAlgorithm, 1)
```
"""
@parameter struct EncodingParameters{C <: Union{AbstractArray{Complex}, Nothing}} <: AbstractMRIRecoEncodingParameters
  correctionMap::C = nothing
  
  method::Union{String, Nothing} = nothing
  toeplitz::Union{Bool, Nothing} = nothing
  oversamplingFactor::Union{Float64, Nothing} = nothing
  kernelSize::Union{Int, Nothing} = nothing
  K::Union{Int, Nothing} = nothing
  K_tol::Union{Float64, Nothing} = nothing
end

"""
    CustomEncodingParameters{E} <: AbstractMRIRecoEncodingParameters

Parameters with user-provided encoding operators.

# Type Parameters
- `E` - The type of the encoding operators

# Fields
- `encodingOps::E` - Pre-computed encoding operators

# Callable Interface
    (ep::CustomEncodingParameters)(::Type{<:AbstractMRIRecoAlgorithm}, slice::Int)

# Returns
- The encoding operators for the given slice

# Example
```julia
ep = CustomEncodingParameters(; encodingOps=myEncodingOps)
encOps = ep(AbstractIterativeMRIRecoAlgorithm, 1)
```
"""
@parameter struct CustomEncodingParameters{E} <: AbstractMRIRecoEncodingParameters
  encodingOps::E
end

function buildEncodingParams(ep::EncodingParameters, S)
  kwargs = (; S = S, copyOpsFn = copyOpsFn(S))
  
  !isnothing(ep.method) && (kwargs = merge(kwargs, (; method = ep.method)))
  !isnothing(ep.toeplitz) && (kwargs = merge(kwargs, (; toeplitz = ep.toeplitz)))
  !isnothing(ep.oversamplingFactor) && (kwargs = merge(kwargs, (; oversamplingFactor = ep.oversamplingFactor)))
  !isnothing(ep.kernelSize) && (kwargs = merge(kwargs, (; kernelSize = ep.kernelSize)))
  !isnothing(ep.K) && (kwargs = merge(kwargs, (; K = ep.K)))
  !isnothing(ep.K_tol) && (kwargs = merge(kwargs, (; K_tol = ep.K_tol)))
  !isnothing(ep.correctionMap) && (kwargs = merge(kwargs, (; correctionMap = ep.correctionMap)))
  
  return kwargs
end

function (ep::EncodingParameters)(::Type{<:AbstractMRIRecoAlgorithm}, slice::Int)
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  
  encParams = buildEncodingParams(ep, S)
  return encodingOps_simple(acqData, reconSize; slice=slice, encParams...)
end

function (ep::CustomEncodingParameters)(::Type{<:AbstractMRIRecoAlgorithm}, slice::Int)
  return ep.encodingOps[:, slice]
end

"""
    SubspaceEncodingParameters{E, B} <: AbstractMRIRecoEncodingParameters

Encoding parameters with subspace projection for dimensionality reduction.

# Type Parameters
- `E` - The inner encoding parameters type
- `B` - The basis matrix type

# Fields
- `inner::E` - Inner encoding parameters
- `basis::B` - Subspace basis matrix

# Callable Interface
    (ep::SubspaceEncodingParameters)(::Type{<:AbstractMRIRecoAlgorithm}, slice::Int, args...)

# Returns
- Encoding operator composed with SubspaceOp

# Example
```julia
ep = SubspaceEncodingParameters(; inner=EncodingParameters(), basis=basisMatrix)
encOp = ep(MultiCoilMultiEchoSubspaceReconstruction, 1, senseMaps)
```
"""
@parameter struct SubspaceEncodingParameters{
    E <: AbstractMRIRecoEncodingParameters,
    B <: AbstractMatrix
} <: AbstractMRIRecoEncodingParameters
  inner::E
  basis::B
end

function (ep::SubspaceEncodingParameters)(algoT::Type{<:AbstractMRIRecoAlgorithm}, slice::Int, args...)
  acqData = ctx_acqData()
  numContr = numContrasts(acqData)
  
  E = ep.inner(algoT, slice, args...)
  subOp = SubspaceOp(ep.basis, acqData.encodingSize, numContr)
  return ∘(E, subOp)
end