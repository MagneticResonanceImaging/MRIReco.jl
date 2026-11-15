export AbstractMRIRecoEncodingParameters
export EncodingParameters, CustomEncodingParameters, SubspaceEncodingParameters

"""
    AbstractMRIRecoEncodingParameters <: AbstractMRIRecoParameters

Abstract base type for encoding operator parameters used in MRI reconstruction.
"""
abstract type AbstractMRIRecoEncodingParameters <: AbstractMRIRecoParameters end

"""
    EncodingParameters{C} <: AbstractMRIRecoEncodingParameters

Parameters for computing encoding operators from acquisition data.

# Type Parameters
- `C` - The type of the correctionMap (Union{AbstractArray{Complex}, Nothing})

# Fields
All fields are `nothing` by default, and the operator then uses its own default given below.
- `correctionMap::C` - Complex correction map `z` for off-resonance and relaxation, in the precision of the data, see
  [Correct Off-Resonance](@ref correct-off-resonance) (default: none)
- `method::Union{String, Nothing}` - Approximation of the term `exp(-z t)` with the sampling times `t`: `"nfft"`
  (default), `"leastsquare"`, `"hist"` or `"svd"`, see [`FieldmapNFFTOp`](@ref)
- `toeplitz::Union{Bool, Nothing}` - Whether the normal operator uses the Toeplitz approach (default: false)
- `oversamplingFactor::Union{Float64, Nothing}` - NFFT oversampling factor (default: 1.25)
- `kernelSize::Union{Int, Nothing}` - NFFT kernel size (default: 3)
- `K::Union{Int, Nothing}` - Number of terms of `"leastsquare"` and `"hist"` and maximum number of terms of `"svd"`
  (default: 20); `"nfft"` determines the number itself
- `K_tol::Union{Float64, Nothing}` - Tolerance with which `"svd"` chooses the number of terms (default: 1e-3)

`method`, `K` and `K_tol` only apply with a correction map; `toeplitz`, `oversamplingFactor` and `kernelSize` only
change non-Cartesian data without a correction map.

# Callable Interface
    (ep::EncodingParameters)(algo, slice::Int)
    (ep::EncodingParameters)(algo, slice::Int, senseMaps)

# Returns
- The encoding operators of the given slice: one per contrast for the reconstructions that treat the contrasts
  separately, one for all contrasts for the others. The second form, used by the reconstructions with coil
  sensitivities, includes the sensitivities.

# Note
- Gets the data, the reconstruction size and the storage type from the [`MRIRecoContext`](@ref) of the reconstruction

# Example
```julia
ep = EncodingParameters()

# With a correction map
ep = EncodingParameters(; correctionMap = cmap, method = "leastsquare", K = 10)
```
"""
@parameter struct EncodingParameters{C <: Union{AbstractArray{<:Complex}, Nothing}} <: AbstractMRIRecoEncodingParameters
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
- `encodingOps::E` - Pre-computed encoding operators of `"standard"` and `"multiCoil"`, a matrix with one row per
  contrast and one column per slice; for `"multiCoil"`, they include the coil sensitivities, as those of
  [`encodingOps_parallel`](@ref). `"multiCoil"` still needs its [`CoilParameters`](@ref), and with `noiseData`, it
  decorrelates the data but not these operators.

# Callable Interface
    (ep::CustomEncodingParameters)(algo, slice::Int, args...)

# Returns
- The encoding operators for the given slice, `encodingOps[:, slice]`

# Example
```julia
ep = CustomEncodingParameters(; encodingOps=myEncodingOps)
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

function (ep::CustomEncodingParameters)(::Type{<:AbstractMRIRecoAlgorithm}, slice::Int, args...)
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
- `basis::B` - Subspace basis matrix with one row per contrast and one column per basis function

# Callable Interface
    (ep::SubspaceEncodingParameters)(algo, slice::Int, args...)

# Returns
- Encoding operator of `inner` composed with a [`SubspaceOp`](@ref) of the basis

# Example
```julia
ep = SubspaceEncodingParameters(; inner=EncodingParameters(), basis=basisMatrix)
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