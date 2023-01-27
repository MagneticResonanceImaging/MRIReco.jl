export defaultRecoParams

function defaultRecoParams()
  params = Dict{Symbol,Any}()
  params[:reco] = "direct"
  params[:sparseTrafoName] = "Wavelet"
  params[:regularization] = "L1"
  params[:λ] = 0.0
  params[:normalizeReg] = false
  params[:solver] = "admm"
  params[:ρ] = 5.e-2
  params[:iterations] = 30

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
# * `noiseData::Array{ComplexF64}`        - noise Acquisition used for noise uncorrelation (pre-whitening)
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
* `reg::Regularization`                 - Regularization to be used
* `normalize::Bool`                     - adjust regularization parameter according to the size of k-space data
* `solvername::String`                  - name of the solver to use
* `senseMaps::Array{Complex{<:AbstractFloat}}`        - coil sensitivities
* `correctionMap::Array{Complex{<:AbstractFloat}}`    - fieldmap for the correction of off-resonance effects
* `method::String="nfft"`               - method to use for time-segmentation when correctio field inhomogeneities

`sparseTrafo` and `reg` can also be speficied using their names in form of a string.
"""
function setupIterativeReco(acqData::AcquisitionData{T}, recoParams::Dict) where T

  red3d = ndims(trajectory(acqData,1))==2 && length(recoParams[:reconSize])==3
  if red3d  # acqData is 3d data converted to 2d
    reconSize = (recoParams[:reconSize][2], recoParams[:reconSize][3])
    recoParams[:reconSize] = reconSize
  else  # regular case
    reconSize = recoParams[:reconSize]
  end

  # density weights
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    weights = samplingDensity(acqData,reconSize)
  else
    numContr = numContrasts(acqData)
    weights = Array{Vector{Complex{T}}}(undef,numContr)
    for contr=1:numContr
      numNodes = size(acqData.kdata[contr],1)
      weights[contr] = [1.0/sqrt(prod(reconSize)) for node=1:numNodes]
    end
  end

  # sparsifying transform
  if haskey(recoParams,:sparseTrafo) && typeof(recoParams[:sparseTrafo]) != String
    sparseTrafo = recoParams[:sparseTrafo]
  elseif haskey(recoParams,:sparseTrafo)
    sparseTrafoName = get(recoParams, :sparseTrafo, "nothing")
    sparseTrafo = SparseOp(Complex{T},sparseTrafoName, reconSize; recoParams...)
  else
    sparseTrafo=nothing
  end

  # bare regularization (without sparsifying transform)
  regName = get(recoParams, :regularization, "L1")
  λ = T(get(recoParams,:λ,0.0))
  reg = Regularization(regName, λ; shape=reconSize, recoParams...)

  # normalize regularizer ?
  normalize = get(recoParams, :normalizeReg, false)

  encOps = get(recoParams, :encodingOps, nothing)

  # solvername
  solvername = get(recoParams, :solver, "fista")

  # sensitivity maps
  senseMaps = get(recoParams, :senseMaps, Complex{T}[])
  if red3d && !isempty(senseMaps) # make sure the dimensions match the trajectory dimensions
    senseMaps = permutedims(senseMaps,[2,3,1,4])
  end

  # noise data acquisition [samples, coils]
  noiseData = get(recoParams, :noiseData, Complex{T}[])
  if isempty(noiseData)
    L_inv = []
  else
    L = cholesky(covariance(noiseData), check = true)
    L_inv = inv(L.L) #noise decorrelation matrix
  end

  return reconSize, weights, L_inv, sparseTrafo, vec(reg), normalize, encOps, solvername, senseMaps
end


function getEncodingOperatorParams(; kargs...)
  encKeys = [:correctionMap, :method, :toeplitz, :oversamplingFactor, :kernelSize, :K, :K_tol]
  return Dict([key=>kargs[key] for key in intersect(keys(kargs),encKeys)])
end

# convenience methods
volumeSize(reconSize::NTuple{2,Int}, numSlice::Int) = (reconSize..., numSlice)
volumeSize(reconSize::NTuple{3,Int}, numSlice::Int) = reconSize
