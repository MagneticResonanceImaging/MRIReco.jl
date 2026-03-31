export AbstractCoilParameters
export CoilParameters

abstract type AbstractCoilParameters <: AbstractMRIRecoParameters end

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

function (params::CoilParameters)(::Type{<:AbstractIterativeMRIRecoAlgorithm})
  acqData = ctx_acqData()
  reconSize = ctx_reconSize()
  arrayType = ctx_arrayType()

  red3d = ndims(trajectory(acqData, 1)) == 2 && length(reconSize) == 3

  senseMaps = params.senseMaps
  if red3d && !isempty(senseMaps)
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