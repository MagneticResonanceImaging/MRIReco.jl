export AbstractDirectMRIRecoParameters

"""
    AbstractDirectMRIRecoParameters <: AbstractMRIRecoParameters

Abstract type for direct (non-iterative) MRI reconstruction parameters.
Performs a direct Fourier-based image reconstruction.
"""
abstract type AbstractDirectMRIRecoParameters <: AbstractMRIRecoParameters end

export DirectMRIReco

"""
    DirectMRIReco{P} <: AbstractDirectMRIRecoAlgorithm

The reconstruction type `"direct"`: a direct image reconstruction which applies the adjoint of the encoding operator
to the density-compensated k-space data.

Its parameter is a [`DirectMRIRecoParameter`](@ref).
"""
@reconstruction mutable struct DirectMRIReco{P <: AbstractDirectMRIRecoParameters} <: AbstractDirectMRIRecoAlgorithm
  @parameter parameter::P
end

export DirectMRIRecoParameter

"""
    DirectMRIRecoParameter{S, C, W, arrT} <: AbstractDirectMRIRecoParameters

Parameters of the direct reconstruction [`DirectMRIReco`](@ref).

# Fields
- `reconSize::S` - Size of the image to reconstruct (default: `nothing`, the encoding size)
- `cmap::C` - Complex correction map for off-resonance and relaxation in the precision of the data, always
  approximated with `"nfft"` (default: `nothing`); in the parameter dictionary of [`reconstruction`](@ref), its key
  is `:cmap`, not `:correctionMap` ([#296](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/296))
- `weights::W` - Has to stay `nothing`: the reconstruction then weights the data with the density compensation
  factors, the squares of the weights from [`samplingDensity`](@ref)
- `arrayType::Type{arrT}` - Array type of the reconstruction (default: `Array`)
"""
@parameter struct DirectMRIRecoParameter{
    S <: Union{Nothing, NTuple{D, Int64} where D},   # User supplied or default recon size
    C <: Union{Nothing, <:AbstractArray{<:Complex}}, # User supplied or no correction map
    W <: Union{Nothing, <:AbstractArray{<:Complex}}, # User supplied or sampling density
    arrT                                             # Array type
  } <: AbstractDirectMRIRecoParameters
  reconSize::S = nothing
  cmap::C = nothing
  weights::W = nothing
  arrayType::Type{arrT} = Array
end


# Setup recon size
(params::DirectMRIRecoParameter{Nothing})(algo::AbstractDirectMRIRecoAlgorithm, acqData::AcquisitionData) = params(algo, acqData, encodingSize(acqData))
(params::DirectMRIRecoParameter)(algo::AbstractDirectMRIRecoAlgorithm, acqData::AcquisitionData) = params(algo, acqData, params.reconSize)

# Prepare image
function (params::AbstractDirectMRIRecoParameters)(algo::AbstractDirectMRIRecoAlgorithm, acqData::AcquisitionData{T}, reconSize::NTuple{D, Int64}) where {T, D}
  encDims = ndims(trajectory(acqData))
  if encDims!=D
    error("reco-dimensionality $D and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)
  Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numChan, numRep)

  Ireco = with(MRIRECO_CONTEXT => MRIRecoContext(reconSize, acqData, params.arrayType)) do
    return params(algo, Ireco, acqData)
  end

  Ireco = reshape(Ireco, volumeSize(reconSize, numSl)..., numContr, numChan, numRep)
  return makeAxisArray(Ireco, acqData)
end

# Reconstruction loop
function (params::DirectMRIRecoParameter)(::Type{<:AbstractDirectMRIRecoAlgorithm}, Ireco::Array{Complex{T}, 5}, acqData::AcquisitionData{T}) where T
  numSl = size(Ireco, 2)
  numContr = size(Ireco, 3)
  numChan = size(Ireco, 4)
  numRep = size(Ireco, 5)

  correctionMap = isnothing(params.cmap) ? similar(Ireco, 0) : params.cmap
  reconSize = ctx_reconSize()
  weights = isnothing(params.weights) ? samplingDensity(acqData, reconSize) : params.weights
  arrayType = ctx_arrayType()
  S = ctx_storageType()

  for i = 1:numSl
    F = encodingOps_simple(acqData, reconSize, slice=i, correctionMap=correctionMap, S = S)
    for k = 1:numContr
      for j = 1:numChan
        for l = 1:numRep
            kdata = arrayType(kData(acqData,k,j,i,rep=l)) .* arrayType((weights[k].^2))
            I = adjoint(F[k]) * kdata
            if isCircular( trajectory(acqData, k) )
              circularShutter!(reshape(I, reconSize), 1.0)
            end
            Ireco[:,i,k,j,l] = Array(I)
        end
      end
    end
  end

  return Ireco
end