export reconstruction_simple, reconstruction_multiEcho, reconstruction_multiCoil, reconstruction_multiCoilMultiEcho, reconstruction_lowRank, RecoParameters

"""
Performs iterative image reconstruction independently for the data of all coils,
contrasts and slices

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Regularization`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
* `solvername::String`                  - name of the solver to use
* (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_simple( acqData::AcquisitionData
                              , reconSize::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = sparseTrafo

  # reconstruction
  Ireco = zeros(ComplexF64, prod(reconSize), numSl, numContr, numChan)
  for k = 1:numSl
    F = encodingOps2d_simple(acqData, reconSize, slice=k, correctionMap=correctionMap, method=method)
    for j = 1:numContr
      W = WeightingOp(weights[j])
      for i = 1:numChan
        kdata = kData(acqData,j,i,k).* weights[j]

        reg2 = deepcopy(reg)
        if normalize
          RegularizedLeastSquares.normalize!(reg2, kdata)
        end
        solver = createLinearSolver(solvername, W*F[j]; reg=reg2, params...)

        I = solve(solver, kdata, startVector=get(params,:startVector,ComplexF64[]),
                                 solverInfo=get(params,:solverInfo,nothing))

        if isCircular( trajectory(acqData, j) )
          circularShutter!(reshape(I, reconSize), 1.0)
        end
        Ireco[:,k,j,i] = I
      end
    end
  end

  Ireco = reshape(Ireco, reconSize..., numSl, numContr, numChan)
  return makeAxisArray(Ireco, acqData)
end

"""
Performs a iterative image reconstruction jointly for all contrasts. Different slices and coil images
are reconstructed independently.

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Regularization`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
* `solvername::String`                  - name of the solver to use
* (`correctionMap::Array{ComplexF64})`  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiEcho(acqData::AcquisitionData
                              , reconSize::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = diagOp( repeat([sparseTrafo],numContr) )

  W = WeightingOp( vcat(weights...) )

  # reconstruction
  Ireco = zeros(ComplexF64, prod(reconSize)*numContr, numChan, numSl)
  for i = 1:numSl
    F = encodingOps2d_multiEcho(acqData, reconSize, slice=k, correctionMap=correctionMap, method=method)
    for j = 1:numChan
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

  Ireco = reshape(Ireco, reconSize..., numContr, numChan, numSl)
  return makeAxisArray(permutedims(Ireco,[1,2,5,3,4]), acqData)
end

"""
Performs a SENSE-type iterative image reconstruction. Different slices and contrasts images
are reconstructed independently.

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Regularization`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
* `solvername::String`                  - name of the solver to use
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
* (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiCoil(acqData::AcquisitionData
                              , reconSize::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , senseMaps::Array{ComplexF64}
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = sparseTrafo

  # solve optimization problem
  Ireco = zeros(ComplexF64, prod(reconSize), numSl, numContr, 1)
  for k = 1:numSl
    E = encodingOps2d_parallel(acqData, reconSize, senseMaps, slice=k, correctionMap=correctionMap, method=method)
    for j = 1:numContr
      W = WeightingOp(weights[j],numChan)
      kdata = multiCoilData(acqData, j, k) .* repeat(weights[j], numChan)

      reg2 = deepcopy(reg)
      if normalize
        RegularizedLeastSquares.normalize!(reg2, kdata)
      end

      solver = createLinearSolver(solvername, W*E[j]; reg=reg2, params...)

      I = solve(solver, kdata)

      if isCircular( trajectory(acqData, j) )
        circularShutter!(reshape(I, reconSize), 1.0)
      end
      Ireco[:,k,j] = I
    end
  end

  Ireco = reshape(Ireco, reconSize..., numSl, numContr, 1)
  return makeAxisArray(Ireco, acqData)
end

"""
Performs a SENSE-type iterative image reconstruction which reconstructs all contrasts jointly.
Different slices are reconstructed independently.

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Regularization`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
* `solvername::String`                  - name of the solver to use
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
* (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)             - method to use for time-segmentation when correctio field inhomogeneities
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiCoilMultiEcho(acqData::AcquisitionData
                              , reconSize::NTuple{2,Int64}
                              , reg::Regularization
                              , sparseTrafo::AbstractLinearOperator
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , senseMaps::Array{ComplexF64}
                              , correctionMap::Array{ComplexF64}=ComplexF64[]
                              , method::String="nfft"
                              , normalize::Bool=false
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}())

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)

  # set sparse trafo in reg
  reg.params[:sparseTrafo] = diagOp( repeat([sparseTrafo],numContr) )

  W = WeightingOp( vcat(weights)..., numChan )

  Ireco = zeros(ComplexF64, prod(reconSize), numContr, numSl)
  for i = 1:numSl
    E = encodingOp_2d_multiEcho_parallel(acqData, reconSize, senseMaps, slice=k, correctionMap=correctionMap, method=method)

    kdata = multiCoilMultiEchoData(acqData, i) .* repeat(weights, numChan)

    reg2 = deepcopy(reg)
    if normalize
      RegularizedLeastSquares.normalize!(reg2, acqData.kdata)
    end
    solver = createLinearSolver(solvername, W*E; reg=reg2, params...)

    Ireco[:,:,i] = solve(solver, kdata)
  end

  Ireco = reshape( permutedims(Ireco, [1,3,2]), recoParams[:reconSize]..., numSl, numContr )
end


"""
Auxilary struct that holds parameters relevant for image reconstruction

# Fields
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `reg::Regularization`                 - Regularization to be used
* `normalize::Bool`                     - adjust regularization parameter according to the size of k-space data
* `solvername::String`                  - name of the solver to use
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
* `correctionMap::Array{ComplexF64}`    - fieldmap for the correction of off-resonance effects
* `method::String="nfft"`               - method to use for time-segmentation when correctio field inhomogeneities
"""
mutable struct RecoParameters{N}
  reconSize::NTuple{N,Int64}
  weights::Vector{Vector{ComplexF64}}
  sparseTrafo::AbstractLinearOperator
  reg::Regularization
  normalize::Bool
  solvername::String
  senseMaps::Array{ComplexF64}
  correctionMap::Array{ComplexF64}
  method::String
end


"""
    setupIterativeReco(acqData::AcquisitionData, recoParams::Dict)

builds relevant parameters and operators from the entries in `recoParams`

# relevant parameters
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `reg::Regularization`                 - Regularization to be used
* `normalize::Bool`                     - adjust regularization parameter according to the size of k-space data
* `solvername::String`                  - name of the solver to use
* `senseMaps::Array{ComplexF64}`        - coil sensitivities
* `correctionMap::Array{ComplexF64}`    - fieldmap for the correction of off-resonance effects
* `method::String="nfft"`               - method to use for time-segmentation when correctio field inhomogeneities

`sparseTrafo` and `reg` can also be speficied using their names in form of a string.
"""
function setupIterativeReco(acqData::AcquisitionData, recoParams::Dict)

  red3d = dims(trajectory(acqData,1))==2 && length(recoParams[:reconSize])==3
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
    weights = Array{Vector{ComplexF64}}(undef,numContr)
    for contr=1:numContr
      numNodes = size(acqData.kdata[contr],1)
      weights[contr] = [1.0/sqrt(prod(reconSize)) for node=1:numNodes]
    end
  end

  # sparsifying transform
  if haskey(recoParams,:sparseTrafo) && typeof(recoParams[:sparseTrafo]) != String
    sparseTrafo = recoParams[:sparseTrafo]
  else
    sparseTrafoName = get(recoParams, :sparseTrafo, "nothing")
    sparseTrafo = SparseOp(sparseTrafoName, reconSize; recoParams...)
  end

  # bare regularization (without sparsifying transform)
  regName = get(recoParams, :regularization, "L1")
  λ = get(recoParams,:λ,0.0)
  reg = Regularization(regName, λ; shape=reconSize, recoParams...)

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

  return RecoParameters(reconSize, weights, sparseTrafo, reg, normalize, solvername, senseMaps, cmap, method)
end
