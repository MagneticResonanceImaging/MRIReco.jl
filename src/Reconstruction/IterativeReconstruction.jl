export reconstruction_simple, reconstruction_multiEcho, reconstruction_multiCoil, reconstruction_multiCoilMultiEcho, reconstruction_lowRank, RecoParameters

"""
Performs iterative image reconstruction independently for the data of all coils,
contrasts and slices

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Vector{<:AbstractRegularization}`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{Complex{<:AbstractFloat}}}` - sampling density of the trajectories in acqData
* `solver::Type{<:AbstractLinearSolver}`                  - name of the solver to use
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_simple( acqData::AcquisitionData{T}
                              ; reconSize::NTuple{D,Int64}
                              , reg::Vector{<:AbstractRegularization}
                              , sparseTrafo
                              , weights::Vector{Vector{Complex{T}}}
                              , solver::Type{<:AbstractLinearSolver}
                              , encodingOps=nothing
                              , params...) where {D, T <: AbstractFloat}

  encDims = ndims(trajectory(acqData))
  if encDims!=length(reconSize)
    error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)

  encParams = getEncodingOperatorParams(;params...)

  # set sparse trafo in reg
  temp = []
  for (i,r) in enumerate(reg)
    trafo = sparseTrafo[i]
    if !isnothing(trafo)
      push!(temp, TransformedRegularization(r, trafo))
    else
      push!(temp, r)
    end
  end
  reg = identity.(temp)

  # reconstruction
  Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numChan, numRep)
  #@floop
  for l = 1:numRep, k = 1:numSl
    if encodingOps!=nothing
      F = encodingOps[:,k]
    else
      F = encodingOps_simple(acqData, reconSize, slice=k; encParams...)
    end
    for j = 1:numContr
      W = WeightingOp(Complex{T}; weights=weights[j])
      for i = 1:numChan
        kdata = kData(acqData,j,i,k,rep=l).* weights[j]
        EFull = ∘(W, F[j])#, isWeighting=true)
        EFullᴴEFull = normalOperator(EFull)
        solv = createLinearSolver(solver, EFull; AᴴA=EFullᴴEFull, reg=reg, params...)

        I = solve(solv, kdata, startVector=get(params,:startVector,Complex{T}[]),
                              solverInfo=get(params,:solverInfo,nothing))

        if isCircular( trajectory(acqData, j) )
          circularShutter!(reshape(I, reconSize), 1.0)
        end
        Ireco[:,k,j,i,l] = I
      end
    end
  end
  Ireco = reshape(Ireco, volumeSize(reconSize, numSl)..., numContr, numChan, numRep)

  return makeAxisArray(Ireco, acqData)
end

"""
Performs a iterative image reconstruction jointly for all contrasts. Different slices and coil images
are reconstructed independently.

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Vector{<:AbstractRegularization}`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{Complex{<:AbstractFloat}}}` - sampling density of the trajectories in acqData
* `solver::Type{<:AbstractLinearSolver}`                  - name of the solver to use
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiEcho(acqData::AcquisitionData{T}
                              ; reconSize::NTuple{D,Int64}
                              , reg::Vector{<:AbstractRegularization}
                              , sparseTrafo
                              , weights::Vector{Vector{Complex{T}}}
                              , solver::Type{<:AbstractLinearSolver}
                              , encodingOps=nothing
                              , params...) where {D , T <: AbstractFloat}

  encDims = ndims(trajectory(acqData))
  if encDims!=length(reconSize)
    error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)
  encParams = getEncodingOperatorParams(;params...)

  # set sparse trafo in reg
  temp = []
  for (i,r) in enumerate(reg)
    trafo = sparseTrafo[i]
    if isnothing(trafo)
      trafo = SparseOp(Complex{T},"nothing", reconSize; params...)
    end
    trafo = DiagOp( repeat([trafo],numContr)... )
    push!(temp, TransformedRegularization(r, trafo))
  end
  reg = identity.(temp)

  W = WeightingOp(Complex{T}; weights=vcat(weights...) )

  # reconstruction
  Ireco = zeros(Complex{T}, prod(reconSize)*numContr, numChan, numSl, numRep)
  @floop for l = 1:numRep, i = 1:numSl
    if encodingOps != nothing
      F = encodingOps[i]
    else
      F = encodingOp_multiEcho(acqData, reconSize, slice=i; encParams...)
    end
    for j = 1:numChan
      kdata = multiEchoData(acqData, j, i,rep=l) .* vcat(weights...)
      EFull = ∘(W, F)#, isWeighting=true)
      EFullᴴEFull = normalOperator(EFull)
      solv = createLinearSolver(solver, EFull; AᴴA=EFullᴴEFull, reg=reg, params...)

      Ireco[:,j,i,l] = solve(solv,kdata; params...)
      # TODO circular shutter
    end
  end

  if encDims==2
    # 2d reconstruction
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], numContr, numChan, numSl,numRep)
    Ireco_ = permutedims(Ireco_, [1,2,5,3,4,6])
  else
    # 3d reconstruction
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numContr, numChan,numRep)
  end

  return makeAxisArray(Ireco_, acqData)
end

"""
Performs a SENSE-type iterative image reconstruction. Different slices and contrasts images
are reconstructed independently.

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Vector{<:AbstractRegularization}`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{Complex{<:AbstractFloat}}}` - sampling density of the trajectories in acqData
* `L_inv::Array{Complex{<:AbstractFloat}}`        - noise decorrelation matrix
* `solver::Type{<:AbstractLinearSolver}`                  - name of the solver to use
* `senseMaps::Array{Complex{<:AbstractFloat}}`        - coil sensitivities
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiCoil(acqData::AcquisitionData{T}
                              ; reconSize::NTuple{D,Int64}
                              , reg::Vector{<:AbstractRegularization}
                              , sparseTrafo
                              , weights::Vector{Vector{Complex{T}}}
                              , L_inv::Union{LowerTriangular{Complex{T}, Matrix{Complex{T}}}, Nothing}
                              , solver::Type{<:AbstractLinearSolver}
                              , senseMaps::Array{Complex{T}}
                              , encodingOps=nothing
                              , params...) where {D , T}

  encDims = ndims(trajectory(acqData))
  if encDims!=length(reconSize)
    error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)
  encParams = getEncodingOperatorParams(;params...)

  # noise decorrelation
  senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan)

  # set sparse trafo in reg
  temp = []
  for (i,r) in enumerate(reg)
    trafo = sparseTrafo[i]
    if !isnothing(trafo)
      push!(temp, TransformedRegularization(r, trafo))
    else
      push!(temp, r)
    end
  end
  reg = identity.(temp)


  # solve optimization problem
  Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numRep)
  let reg = reg # Fix @floop warning due to conditional/multiple assignment to reg
    @floop for l = 1:numRep, k = 1:numSl
      if encodingOps != nothing
        E = encodingOps[:,k]
      else
        E = encodingOps_parallel(acqData, reconSize, senseMapsUnCorr; slice=k, encParams...)
      end

      for j = 1:numContr
        W = WeightingOp(Complex{T}; weights=weights[j], rep=numChan)
        kdata = multiCoilData(acqData, j, k, rep=l) .* repeat(weights[j], numChan)
        if !isnothing(L_inv)
          kdata = vec(reshape(kdata, :, numChan) * L_inv')
        end

        EFull = ∘(W, E[j], isWeighting=true)
        EFullᴴEFull = normalOperator(EFull)
        solv = createLinearSolver(solver, EFull; AᴴA=EFullᴴEFull, reg=reg, params...)
        I = solve(solv, kdata; params...)

        if isCircular( trajectory(acqData, j) )
          circularShutter!(reshape(I, reconSize), 1.0)
        end
        Ireco[:,k,j,l] = I
      end
    end
  end
  Ireco_ = reshape(Ireco, volumeSize(reconSize, numSl)..., numContr, 1,numRep)

  return makeAxisArray(Ireco_, acqData)
end


"""
Performs a SENSE-type iterative image reconstruction which reconstructs all contrasts jointly.
Different slices are reconstructed independently.

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Vector{<:AbstractRegularization}`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{Complex{<:AbstractFloat}}}` - sampling density of the trajectories in acqData
* `solver::Type{<:AbstractLinearSolver}`                  - name of the solver to use
* `senseMaps::Array{Complex{<:AbstractFloat}}`        - coil sensitivities
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiCoilMultiEcho(acqData::AcquisitionData{T}
                              ; reconSize::NTuple{D,Int64}
                              , reg::Vector{<:AbstractRegularization}
                              , sparseTrafo
                              , weights::Vector{Vector{Complex{T}}}
                              , L_inv::Union{LowerTriangular{Complex{T}, Matrix{Complex{T}}}, Nothing}
                              , solver::Type{<:AbstractLinearSolver}
                              , senseMaps::Array{Complex{T}}
                              , encodingOps=nothing
                              , params...) where {D, T}

  encDims = ndims(trajectory(acqData))
  if encDims!=length(reconSize)
    error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)
  encParams = getEncodingOperatorParams(;params...)

  # noise decorrelation
  senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan)

  # set sparse trafo in reg
  temp = []
  for (i,r) in enumerate(reg)
    trafo = sparseTrafo[i]
    if isnothing(trafo)
      trafo = SparseOp(Complex{T},"nothing", reconSize; params...)
    end
    trafo = DiagOp( repeat([trafo],numContr)... )
    push!(temp, TransformedRegularization(r, trafo))
  end
  reg = identity.(temp)

  W = WeightingOp(Complex{T}; weights=vcat(weights...), rep=numChan )

  Ireco = zeros(Complex{T}, prod(reconSize)*numContr, numSl, numRep)
  @floop for l = 1:numRep, i = 1:numSl
    if encodingOps != nothing
      E = encodingOps[i]
    else
      E = encodingOp_multiEcho_parallel(acqData, reconSize, senseMapsUnCorr; slice=i, encParams...)
    end

    kdata = multiCoilMultiEchoData(acqData, i) .* repeat(vcat(weights...), numChan)
    if !isnothing(L_inv)
      kdata = vec(reshape(kdata, :, numChan) * L_inv')
    end

    EFull = ∘(W, E)#, isWeighting=true)
    EFullᴴEFull = normalOperator(EFull)
    solv = createLinearSolver(solver, EFull; AᴴA=EFullᴴEFull, reg=reg, params...)

    Ireco[:,i,l] = solve(solv, kdata; params...)
  end


  if encDims==2
    # 2d reconstruction
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], numContr, numSl, 1, numRep)
    Ireco_ = permutedims(Ireco_, [1,2,4,3,5,6])
  else
    # 3d reconstruction
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numContr, 1, numRep)
  end

  return makeAxisArray(Ireco_, acqData)
end



"""
Performs a SENSE-type iterative image reconstruction which reconstructs all contrasts jointly.
Different slices are reconstructed independently.

# Arguments
* `acqData::AcquisitionData`            - AcquisitionData object
* `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
* `reg::Vector{<:AbstractRegularization}`                 - Regularization to be used
* `sparseTrafo::AbstractLinearOperator` - sparsifying transformation
* `weights::Vector{Vector{Complex{<:AbstractFloat}}}` - sampling density of the trajectories in acqData
* `solver::Type{<:AbstractLinearSolver}`                  - name of the solver to use
* `senseMaps::Array{Complex{<:AbstractFloat}}`        - coil sensitivities
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiCoilMultiEcho_subspace(acqData::AcquisitionData{T}
                              ; reconSize::NTuple{D,Int64}
                              , reg::Vector{<:AbstractRegularization}
                              , sparseTrafo
                              , weights::Vector{Vector{Complex{T}}}
                              , solver::Type{<:AbstractLinearSolver}
                              , senseMaps::Array{Complex{T}}
                              , encodingOps=nothing
                              , params...) where {D, T}

  encDims = ndims(trajectory(acqData))
  if encDims!=length(reconSize)
    error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)
  numBasis = size(params[:basis],2)
   
  encParams = getEncodingOperatorParams(;params...)

  # set sparse trafo in reg
  temp = []
  for (i,r) in enumerate(reg)
    trafo = sparseTrafo[i]
    if isnothing(trafo)
      trafo = SparseOp(Complex{T},"nothing", reconSize; params...)
    end
    trafo = DiagOp( repeat([trafo],numContr)... )
    push!(temp, TransformedRegularization(r, trafo))
  end
  reg = identity.(temp)

  W = WeightingOp(Complex{T}; weights=vcat(weights...), rep=numChan )

  Ireco = zeros(Complex{T}, prod(reconSize)*numBasis, numSl, numRep)
  @floop for l = 1:numRep, i = 1:numSl
    if encodingOps != nothing
      E = encodingOps[i]
    else
      E = encodingOp_multiEcho_parallel(acqData, reconSize, senseMaps; slice=i, encParams...)
      subOp = SubspaceOp(params[:basis],acqData.encodingSize,numContr)
      E = ∘(E,subOp)
    end

    kdata = multiCoilMultiEchoData(acqData, i) .* repeat(vcat(weights...), numChan)

    EFull = ∘(W, E)#, isWeighting=true)
    EFullᴴEFull = normalOperator(EFull)
    solv = createLinearSolver(solver, EFull; AᴴA=EFullᴴEFull, reg=reg, params...)

    Ireco[:,i,l] = solve(solv, kdata; params...)
  end


  if encDims==2
    # 2d reconstruction
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], numBasis, numSl, 1, numRep)
    Ireco_ = permutedims(Ireco_, [1,2,4,3,5,6])
  else
    # 3d reconstruction
    Ireco_ = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numBasis, 1, numRep)
  end

  return makeAxisArray(Ireco_, acqData)
end