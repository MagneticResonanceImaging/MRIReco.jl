export reconstruction_simple, reconstruction_multiEcho, reconstruction_multiCoil, reconstruction_multiCoilMultiEcho, reconstruction_lowRank, RecoParameters

"""
  CS-Sense Reconstruction using sparsity in the wavelet domain
"""
function reconstruction_simple( acqData::AcquisitionData
                              , shape::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = sparseTrafo

  # reconstruction
  Ireco = zeros(ComplexF64, prod(shape), acqData.numSlices, acqData.numEchoes, acqData.numCoils)
  for k = 1:acqData.numSlices
    F = encodingOps2d_simple(acqData, shape, slice=k, correctionMap=correctionMap, method=method)
    for j = 1:acqData.numEchoes
      W = WeightingOp(weights[j])
      for i = 1:acqData.numCoils
        kdata = kData(acqData,j,i,k).* weights[j]

        reg2 = deepcopy(reg)
        if normalize
          RegularizedLeastSquares.normalize!(reg2, kdata)
        end
        solver = createLinearSolver(solvername, W*F[j]; reg=reg2, params...)

        I = solve(solver, kdata)

        if isCircular( trajectory(acqData, j) )
          circularShutter!(reshape(I, shape), 1.0)
        end
        Ireco[:,k,j,i] = I
      end
    end
  end

  Ireco = reshape(Ireco, shape..., acqData.numSlices, acqData.numEchoes, acqData.numCoils)
  return makeAxisArray(Ireco, acqData)
end

"""
  CS Reconstruction using a joint encoding operator for the different echos
  and regularization on the multi-echo data
"""
function reconstruction_multiEcho(acqData::AcquisitionData
                              , shape::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = diagOp( repeat([sparseTrafo],acqData.numEchoes) )

  W = WeightingOp( vcat(weights...) )

  # reconstruction
  Ireco = zeros(ComplexF64, prod(shape)*acqData.numEchoes, acqData.numCoils, acqData.numSlices)
  for i = 1:acqData.numSlices
    F = encodingOps2d_multiEcho(acqData, shape, slice=k, correctionMap=correctionMap, method=method)
    for j = 1:acqData.numCoils
      kdata = multiEchoData(acqData, j, i) .* weights

      reg2 = deepcopy(reg)
      if normalize
        RegularizedLeastSquares.normalize!(reg2, kdata)
      end
      solver = createLinearSolver(solvername, W*F; reg=reg2, params...)

      Ireco[:,j,i] = solve(solver,kdata)
      # TODO circular shutter
    end
  end

  Ireco = reshape(Ireco, shape..., acqData.numEchoes, acqData.numCoils, acqData.numSlices)
  return makeAxisArray(permutedims(Ireco,[1,2,5,3,4]), acqData)
end

"""
  CS-Sense Reconstruction
"""
function reconstruction_multiCoil(acqData::AcquisitionData
                              , shape::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , senseMaps::Array{ComplexF64}
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = sparseTrafo

  # solve optimization problem
  Ireco = zeros(ComplexF64, prod(shape), acqData.numSlices, acqData.numEchoes, 1)
  for k = 1:acqData.numSlices
    E = encodingOps2d_parallel(acqData, shape, senseMaps, slice=k, correctionMap=correctionMap, method=method)
    for j = 1:acqData.numEchoes
      W = WeightingOp(weights[j],acqData.numCoils)
      kdata = multiCoilData(acqData, j, k) .* repeat(weights[j], acqData.numCoils)

      reg2 = deepcopy(reg)
      if normalize
        RegularizedLeastSquares.normalize!(reg2, kdata)
      end

      solver = createLinearSolver(solvername, W*E[j]; reg=reg2, params...)

      I = solve(solver, kdata)

      if isCircular( trajectory(acqData, j) )
        circularShutter!(reshape(I, shape), 1.0)
      end
      Ireco[:,k,j] = I
    end
  end

  Ireco = reshape(Ireco, shape..., acqData.numSlices, acqData.numEchoes, 1)
  return makeAxisArray(Ireco, acqData)
end

function reconstruction_multiCoilMultiEcho(acqData::AcquisitionData
                              , shape::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , senseMaps::Array{ComplexF64}
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = diagOp( repeat([sparseTrafo],acqData.numEchoes) )

  W = WeightingOp( vcat(weights)..., acqData.numCoils )

  Ireco = zeros(ComplexF64, prod(shape), acqData.numEchoes, acqData.numSlices)
  for i = 1:acqData.numSlices
    E = encodingOp_2d_multiEcho_parallel(acqData, shape, senseMaps, slice=k, correctionMap=correctionMap, method=method)

    kdata = multiCoilMultiEchoData(acqData, i) .* repeat(weights, acqData.numCoils)

    reg2 = deepcopy(reg)
    if normalize
      RegularizedLeastSquares.normalize!(reg2, acqData.kdata)
    end
    solver = createLinearSolver(solvername, W*E; reg=reg2, params...)

    Ireco[:,:,i] = solve(solver, kdata)
  end

  Ireco = reshape( permutedims(Ireco, [1,3,2]), recoParams[:shape]..., acqData.numSlices, acqData.numEchoes )
end


###########################################################################
# setup regularization, sparsifying transform  and density weights for reco
###########################################################################
mutable struct RecoParameters{N}
  shape::NTuple{N,Int64}
  weights::Vector{Vector{ComplexF64}}
  sparseTrafo::AbstractLinearOperator
  reg::Regularization
  normalize::Bool
  solvername::String
  senseMaps::Array{ComplexF64}
  correctionMap::Array{ComplexF64}
  method::String
end

function setupIterativeReco(acqData::AcquisitionData, recoParams::Dict)

  red3d = dims(trajectory(acqData,1))==2 && length(recoParams[:shape])==3
  if red3d  # acqData is 3d data converted to 2d
    shape = (recoParams[:shape][2], recoParams[:shape][3])
    recoParams[:shape] = shape
  else  # regular case
    shape = recoParams[:shape]
  end

  # density weights
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    weights = samplingDensity(acqData,shape)
  else
    weights = [1.0/sqrt(prod(shape)) for echo=1:acqData.numEcoes]
  end

  # sparsifying transform
  if haskey(recoParams,:sparseTrafo)
    sparseTrafo = recoParams[:sparseTrafo]
  else
    sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
    sparseTrafo = SparseOp(sparseTrafoName, shape; recoParams...)
  end

  # bare regularization (without sparsifying transform)
  regName = get(recoParams, :regularization, "L1")
  λ = get(recoParams,:λ,0.0)
  reg = Regularization(regName, λ; recoParams...)

  # normalize regularizer ?
  normalize = get(recoParams, :normalizeReg, false)

  # solvername
  solvername = get(recoParams, :solver, "fista")

  # sensitivity maps
  senseMaps = get(recoParams, :senseMaps, ComplexF64[])
  if red3d && !isempty(senseMaps) # make sure the dimensions match the trajectory dimensions
    senseMaps = permutedims(senseMaps,[2,3,1,4])
  end

  # field map
  cmap = get(recoParams, :correctionMap, ComplexF64[])
  method = get(recoParams, :method, "nfft")

  return RecoParameters(shape, weights, sparseTrafo, reg, normalize, solvername, senseMaps, cmap, method)
end
