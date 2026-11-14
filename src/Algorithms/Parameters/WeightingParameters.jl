export AbstractMRIRecoWeightingParameters
export DensityWeightingParameters, UniformWeightingParameters, CustomWeightingParameters

"""
    AbstractMRIRecoWeightingParameters <: AbstractMRIRecoParameters

Abstract base type for weighting parameters used in MRI reconstruction.
Subtypes define how weights are computed for the reconstruction.
"""
abstract type AbstractMRIRecoWeightingParameters <: AbstractMRIRecoParameters end

"""
    DensityWeightingParameters <: AbstractMRIRecoWeightingParameters

Weighting parameters that compute density compensation weights from the acquisition data.

# Description
Uses the trajectory information in the acquisition data to compute appropriate
density compensation weights. This is the default weighting scheme for most
reconstructions.

# Callable Interface
    (wp::DensityWeightingParameters)(::Type{<:AbstractMRIRecoAlgorithm})

# Returns
- `Vector{Vector{Complex{T}}}` - One weight vector per contrast

# Example
```julia
wp = DensityWeightingParameters()
weights = wp(AbstractIterativeMRIRecoAlgorithm)
```
"""
@parameter struct DensityWeightingParameters <: AbstractMRIRecoWeightingParameters
end

"""
    UniformWeightingParameters <: AbstractMRIRecoWeightingParameters

Weighting parameters that use uniform weights for all k-space samples.

# Description
Assigns uniform weight `1/sqrt(prod(reconSize))` to each k-space sample.

# Callable Interface
    (wp::UniformWeightingParameters)(::Type{<:AbstractMRIRecoAlgorithm})

# Returns
- `Vector{Vector{Complex{T}}}` - One weight vector per contrast

# Example
```julia
wp = UniformWeightingParameters()
weights = wp(AbstractIterativeMRIRecoAlgorithm)
```
"""
@parameter struct UniformWeightingParameters <: AbstractMRIRecoWeightingParameters  
end

"""
    CustomWeightingParameters{W} <: AbstractMRIRecoWeightingParameters

Weighting parameters that use user-provided custom weights.

# Type Parameters
- `W` - The type of the weights (typically `Vector{Vector{Complex{<:AbstractFloat}}}`)

# Fields
- `weights::W` - User-provided weights

# Description
Allows explicit specification of weighting vectors. The weights should be
provided as a vector of vectors, one per contrast.

# Callable Interface
    (wp::CustomWeightingParameters)(::Type{<:AbstractMRIRecoAlgorithm})

# Returns
- The stored `weights` vector

# Validation
Weights are validated during reconstruction to ensure they match the
acquisition data structure (same number of contrasts and samples).

# Example
```julia
weights = [[1.0+0im, 2.0+0im], [3.0+0im, 4.0+0im]]
wp = CustomWeightingParameters(; weights=weights)
result = wp(AbstractIterativeMRIRecoAlgorithm)
```
"""
@parameter struct CustomWeightingParameters{W} <: AbstractMRIRecoWeightingParameters
  weights::W
end

function (wp::DensityWeightingParameters)(::Type{<:AbstractMRIRecoAlgorithm})
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  return map(weight -> S(weight), samplingDensity(acqData, reconSize))
end

function (wp::UniformWeightingParameters)(::Type{<:AbstractMRIRecoAlgorithm})
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  T = real(eltype(S))
  
  numContr = numContrasts(acqData)
  weights = Vector{eltype(S)}[]
  for contr in 1:numContr
    numNodes = length(acqData.kdata[contr])
    push!(weights, S([one(T)/sqrt(prod(reconSize)) for _ in 1:numNodes]))
  end
  return weights
end

function (wp::CustomWeightingParameters)(::Type{<:AbstractMRIRecoAlgorithm})
  acqData = ctx_acqData()
  
  # Validation: check number of contrasts
  numContr_expected = numContrasts(acqData)
  numContr_provided = length(wp.weights)
  if numContr_provided != numContr_expected
    throw(DimensionMismatch("CustomWeightingParameters: expected $numContr_expected contrast weight vectors, got $numContr_provided"))
  end
  
  # Validation: check number of samples per contrast
  for (i, weight_vec) in enumerate(wp.weights)
    numSamples_expected = length(acqData.kdata[i])
    numSamples_provided = length(weight_vec)
    if numSamples_provided != numSamples_expected
      throw(DimensionMismatch("CustomWeightingParameters: contrast $i expected $numSamples_expected samples, got $numSamples_provided"))
    end
  end
  
  return wp.weights
end