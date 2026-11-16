export AbstractCoilParameters
export CoilParameters

abstract type AbstractCoilParameters <: AbstractMRIRecoParameters end

"""
    CoilParameters{T, S, N} <: AbstractCoilParameters
    CoilParameters(; senseMaps = nothing, noiseData = nothing)

Coil sensitivities and noise data of the SENSE-type reconstructions `"multiCoil"`, `"multiCoilMultiEcho"` and
`"multiCoilMultiEchoSubspace"`, as complex arrays in the precision of the data. The constructor needs at least one
of the two, and the reconstructions with [`EncodingParameters`](@ref) need `senseMaps`.

# Fields
* `senseMaps::S` - Coil sensitivities with the dimensions x, y, z or slice, and coil, as `espirit` returns them when
                   called on an `AcquisitionData` with `nmaps = 1`, see
                   [Estimate Coil Sensitivities with ESPIRiT](@ref estimate-coil-sensitivities)
* `noiseData::N` - Noise samples as a matrix of samples × coils; if given, the k-space data and the coil
                   sensitivities are decorrelated, see [Noise Decorrelation](@ref)
"""
@parameter constructor = false struct CoilParameters{
    T,
    S <: AbstractArray{Complex{T}},
    N <: AbstractArray{Complex{T}}
  } <: AbstractCoilParameters
  senseMaps::S
  noiseData::N
end
CoilParameters(; senseMaps = nothing, noiseData = nothing) = CoilParameters(senseMaps, noiseData)
CoilParameters(senseMaps::AbstractArray{T}, noiseData::Nothing) where T = CoilParameters(senseMaps, T[])
CoilParameters(senseMaps::Nothing, noiseData::AbstractArray{T}) where T = CoilParameters(T[], noiseData)
CoilParameters(::Nothing, ::Nothing) = throw(ArgumentError("Cannot construct CoilParameters with no senseMaps and no noiseData"))

function (params::CoilParameters)(::Type{<:AbstractIterativeMRIRecoAlgorithm})
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()

  senseMaps = arrayType(params.senseMaps)
  if ctx_red3d() && !isempty(senseMaps)
    senseMaps = permutedims(senseMaps, [2, 3, 1, 4])
  end

  numChan = isempty(senseMaps) ? 0 : size(senseMaps, ndims(senseMaps))

  noiseData = arrayType(params.noiseData)
  if isempty(noiseData)
    L_inv = nothing
  else
    psi = convert(typeof(noiseData), covariance(noiseData))
    # Check if approx. hermitian and avoid check in cholesky
    # GPU arrays can have float rounding differences
    @assert isapprox(adjoint(psi), psi)
    L = cholesky(psi; check = false)
    L_inv = inv(L.L) #noise decorrelation matrix
  end

  if !isnothing(L_inv) && !isempty(senseMaps)
    senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan)
  else
    senseMapsUnCorr = senseMaps
  end

  return L_inv, senseMapsUnCorr
end