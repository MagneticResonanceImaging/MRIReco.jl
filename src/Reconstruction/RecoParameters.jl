export defaultRecoParams

function defaultRecoParams()
  params = Dict{Symbol,Any}()
  params[:reco] = "direct"
  params[:sparseTrafoName] = "Wavelet"
  params[:reg] = L1Regularization(0.0)
  params[:normalizeReg] = NoNormalization()
  params[:solver] = ADMM
  params[:ρ] = 5.e-2
  params[:iterations] = 30
  params[:arrayType] = Array

  return params
end

function setupDirectReco(acqData::AcquisitionData{T}, recoParams::Dict) where T
  reconSize = recoParams[:reconSize]
  weights = samplingDensity(acqData, recoParams[:reconSize])
  # field map
  cmap = get(recoParams, :cmap, Complex{T}[])

  return reconSize, weights, cmap
end


# """
# Auxilary struct that holds parameters relevant for image reconstruction

# # Fields
# * `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
# * `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
# * `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
# * `reg::Regularization`                 - Regularization to be used
# * `normalize::Bool`                     - adjust regularization parameter according to the size of k-space data
# * `solvername::String`                  - name of the solver to use
# * `senseMaps::Array{ComplexF64}`        - coil sensitivities
# * `correctionMap::Array{ComplexF64}`    - fieldmap for the correction of off-resonance effects
# * `method::String="nfft"`               - method to use for time-segmentation when correctio field inhomogeneities
# * `noiseData::Array{ComplexF64}`        - noise acquisition for noise decorelation
# """
# mutable struct RecoParameters{N}
#   reconSize::NTuple{N,Int64}
#   weights::Vector{Vector{ComplexF64}}
#   sparseTrafo::AbstractLinearOperator
#   reg::Vector{Regularization}
#   normalize::Bool
#   encodingOps::
#   solvername::String
#   senseMaps::Array{ComplexF64}
#   noiseData::Array{ComplexF64}
# end


"""
    setupIterativeReco(acqData::AcquisitionData, recoParams::Dict)

builds relevant parameters and operators from the entries in `recoParams`

# relevant parameters
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `weights::Vector{Vector{Complex{<:AbstractFloat}}}` - sampling density of the trajectories in acqData
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `reg::AbstractRegularization`                 - Regularization to be used
* `normalize::Bool`                     - adjust regularization parameter according to the size of k-space data
* `solver::Type{<:AbstractLinearSolver}`                  - solver to use
* `senseMaps::Array{Complex{<:AbstractFloat}}`        - coil sensitivities
* `correctionMap::Array{Complex{<:AbstractFloat}}`    - fieldmap for the correction of off-resonance effects
* `method::String="nfft"`               - method to use for time-segmentation when correctio field inhomogeneities
* `noiseData::Array{ComplexF64}`        - noise acquisition for noise decorelation

`sparseTrafo` and `reg` can also be speficied using their names in form of a string.
"""
function setupIterativeReco!(acqData::AcquisitionData{T}, recoParams::Dict) where T

  red3d = ndims(trajectory(acqData,1))==2 && length(recoParams[:reconSize])==3
  if red3d  # acqData is 3d data converted to 2d
    reconSize = (recoParams[:reconSize][2], recoParams[:reconSize][3])
    recoParams[:reconSize] = reconSize
  else  # regular case
    reconSize = recoParams[:reconSize]
  end

  arrayType = get(recoParams, :arrayType, Array)
  vecTc = arrayType{Complex{T}, 1}

  # density weights
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    weights = map(weight -> vecTc(weight), samplingDensity(acqData,reconSize))
  else
    numContr = numContrasts(acqData)
    weights = Array{vecTc}(undef,numContr)
    for contr=1:numContr
      numNodes = size(acqData.kdata[contr],1)
      weights[contr] = vecTc([1.0/sqrt(prod(reconSize)) for node=1:numNodes])
    end
  end

  # bare regularization (without sparsifying transform)
  reg = get(recoParams,:reg,L1Regularization(zero(T)))
  reg = vec(reg)

  # sparsifying transform
  if haskey(recoParams, :sparseTrafo)
    sparseTrafos = recoParams[:sparseTrafo]
    if !(typeof(sparseTrafos) <: Vector)
      sparseTrafos = [sparseTrafos]
    end

    # Construct SparseOp for each string instance
    sparseTrafos = map(x-> x isa String ? SparseOp(Complex{T}, x, reconSize; S = vecTc, recoParams...) : x, sparseTrafos)

    # Fill up SparseOps for remaining reg terms with nothing
    temp = Union{Nothing, eltype(sparseTrafos)}[nothing for i = 1:length(reg)]
    for (i,sparseTrafo) in enumerate(sparseTrafos)
      temp[i] = sparseTrafo
    end

    sparseTrafo = identity.(temp)
  else
    sparseTrafo = fill(nothing, length(reg))
  end

  # normalize regularizer ?
  normalize = get(recoParams, :normalizeReg, NoNormalization())

  encOps = get(recoParams, :encodingOps, nothing)

  # solvername
  solver = get(recoParams, :solver, FISTA)

  # sensitivity maps
  senseMaps = get(recoParams, :senseMaps, vecTc())
  if red3d && !isempty(senseMaps) # make sure the dimensions match the trajectory dimensions
    senseMaps = permutedims(senseMaps,[2,3,1,4])
  end

  # noise data acquisition [samples, coils]
  noiseData = get(recoParams, :noiseData, vecTc())
  if isempty(noiseData)
    L_inv = nothing
  else
    psi = convert(typeof(noiseData), covariance(noiseData))
    L = cholesky(psi, check = true)
    L_inv = inv(L.L) #noise decorrelation matrix
  end


  recoParams[:reconSize] = reconSize
  recoParams[:weights] = weights
  recoParams[:L_inv] = L_inv
  recoParams[:sparseTrafo] = sparseTrafo
  recoParams[:reg] = reg
  recoParams[:normalize] = normalize 
  recoParams[:encOps] = encOps
  recoParams[:solver] = solver
  recoParams[:senseMaps] = senseMaps
  return recoParams
end


function getEncodingOperatorParams(; kargs...)
  encKeys = [:correctionMap, :method, :toeplitz, :oversamplingFactor, :kernelSize, :K, :K_tol]
  return Dict([key=>kargs[key] for key in intersect(keys(kargs),encKeys)])
end

# convenience methods
volumeSize(reconSize::NTuple{2,Int}, numSlice::Int) = (reconSize..., numSlice)
volumeSize(reconSize::NTuple{3,Int}, numSlice::Int) = reconSize

executor(::Type{<:AbstractArray}) = nothing