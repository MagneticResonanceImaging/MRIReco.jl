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
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_simple( acqData::AcquisitionData
                              , reconSize::NTuple{D,Int64}
                              , reg::Vector{Regularization}
                              , sparseTrafo::Trafo
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , normalize::Bool=false
                              , encodingOps=nothing
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}()) where D

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)

  encParams = getEncodingOperatorParams(;params...)

  # set sparse trafo in reg
  reg[1].params[:sparseTrafo] = sparseTrafo

  # reconstruction
  Ireco = zeros(ComplexF64, prod(reconSize), numSl, numContr, numChan)
  @sync for k = 1:numSl
    Threads.@spawn begin
      if encodingOps!=nothing
        F = encodingOps[:,k]
      else
        F = encodingOps_simple(acqData, reconSize, slice=k; encParams...)
      end
      for j = 1:numContr
        W = WeightingOp(weights[j])
        for i = 1:numChan
          kdata = kData(acqData,j,i,k).* weights[j]
          solver = createLinearSolver(solvername, W∘F[j]; reg=reg, params...)

          I = solve(solver, kdata, startVector=get(params,:startVector,ComplexF64[]),
                                 solverInfo=get(params,:solverInfo,nothing))

          if isCircular( trajectory(acqData, j) )
            circularShutter!(reshape(I, reconSize), 1.0)
          end
          Ireco[:,k,j,i] = I
        end
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
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiEcho(acqData::AcquisitionData
                              , reconSize::NTuple{D,Int64}
                              , reg::Vector{Regularization}
                              , sparseTrafo::Trafo
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , normalize::Bool=false
                              , encodingOps=nothing
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}()) where D

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)
  encParams = getEncodingOperatorParams(;params...)

  # set sparse trafo in reg
  reg[1].params[:sparseTrafo] = diagOp( repeat([sparseTrafo],numContr)... )

  W = WeightingOp( vcat(weights...) )

  # reconstruction
  Ireco = zeros(ComplexF64, prod(reconSize)*numContr, numChan, numSl)
  @sync for i = 1:numSl
    Threads.@spawn begin
      if encodingOps != nothing
        F = encodingOps[i]
      else
        F = encodingOp_multiEcho(acqData, reconSize, slice=i; encParams...)
      end
      for j = 1:numChan
        kdata = multiEchoData(acqData, j, i) .* vcat(weights...)
        solver = createLinearSolver(solvername, W∘F; reg=reg, params...)

        Ireco[:,j,i] = solve(solver,kdata; params...)
        # TODO circular shutter
      end
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
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiCoil(acqData::AcquisitionData
                              , reconSize::NTuple{D,Int64}
                              , reg::Vector{Regularization}
                              , sparseTrafo::Trafo
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , senseMaps::Array{ComplexF64}
                              , normalize::Bool=false
                              , encodingOps=nothing
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}()) where D

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)
  encParams = getEncodingOperatorParams(;params...)

  # set sparse trafo in reg
  reg[1].params[:sparseTrafo] = sparseTrafo

  # solve optimization problem
  Ireco = zeros(ComplexF64, prod(reconSize), numSl, numContr, 1)
  @sync for k = 1:numSl
    Threads.@spawn begin
      if encodingOps != nothing
        E = encodingOps[:,k]
      else
        E = encodingOps_parallel(acqData, reconSize, senseMaps; slice=k, encParams...)
      end

      for j = 1:numContr
        W = WeightingOp(weights[j],numChan)
        kdata = multiCoilData(acqData, j, k) .* repeat(weights[j], numChan)

        EFull = ∘(W, E[j], isWeighting=true)
        solver = createLinearSolver(solvername, EFull; reg=reg, params...)
        I = solve(solver, kdata; params...)

        if isCircular( trajectory(acqData, j) )
          circularShutter!(reshape(I, reconSize), 1.0)
        end
        Ireco[:,k,j] = I
      end
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
* (`normalize::Bool=false`)             - adjust regularization parameter according to the size of k-space data
* (`params::Dict{Symbol,Any}`)          - Dict with additional parameters
"""
function reconstruction_multiCoilMultiEcho(acqData::AcquisitionData
                              , reconSize::NTuple{D,Int64}
                              , reg::Vector{Regularization}
                              , sparseTrafo::Trafo
                              , weights::Vector{Vector{ComplexF64}}
                              , solvername::String
                              , senseMaps::Array{ComplexF64}
                              , normalize::Bool=false
                              , encodingOps=nothing
                              , params::Dict{Symbol,Any}=Dict{Symbol,Any}()) where D

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)
  encParams = getEncodingOperatorParams(;params...)

  # set sparse trafo in reg
  reg[1].params[:sparseTrafo] = diagOp( repeat([sparseTrafo],numContr)... )

  W = WeightingOp( vcat(weights...), numChan )

  Ireco = zeros(ComplexF64, prod(reconSize), numContr, numSl)
  @sync for i = 1:numSl
    Threads.@spawn begin
      if encodingOps != nothing
        E = encodingOps[i]
      else
        E = encodingOp_multiEcho_parallel(acqData, reconSize, senseMaps; slice=i, encParams...)
      end

      kdata = multiCoilMultiEchoData(acqData, i) .* repeat(vcat(weights...), numChan)
      solver = createLinearSolver(solvername, W∘E; reg=reg, params...)

      Ireco[:,:,i] = solve(solver, kdata; params...)
    end
  end
  return reshape( permutedims(Ireco, [1,3,2]), reconSize..., numSl, numContr )
end

