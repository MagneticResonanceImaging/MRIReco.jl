export reconstruction_simple, reconstruction_multiEcho, reconstruction_multiCoil, reconstruction_multiCoilMultiEcho, reconstruction_lowRank

"""
  CS-Sense Reconstruction using sparsity in the wavelet domain
"""
function reconstruction_simple(acqData::AcquisitionData, recoParams::Dict)

  # sparsifying transform
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # weight kdata with sampling density
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    weights = samplingDensity(acqData,recoParams[:shape])
  else
    weights = [1.0/sqrt(prod(recoParams[:shape])) for echo=1:acqData.numEcoes]
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  λ = get(recoParams,:λ,0.0)
  normalize = get(recoParams, :normalizeReg, false)

  # reconstruction
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), acqData.numSlices, acqData.numEchoes, acqData.numCoils)
  solvername = get(recoParams, :solver, "fista")

  for k = 1:acqData.numSlices
    F = encodingOps2d_simple(acqData, recoParams, slice=k)
    for j = 1:acqData.numEchoes
      W = WeightingOp(weights[j])
      for i = 1:acqData.numCoils
        reg = Regularization(regName, λ; recoParams...)
        if normalize
          RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
        end
        solver = createLinearSolver(solvername, W*F[j]; reg=reg, recoParams...)
        kdata = kData(acqData,j,i,k).* weights[j] #kData(acqData2,j,i,k)
        I = solve(solver, kdata)

        if isCircular( trajectory(acqData, j) )
          circularShutter!(reshape(I, recoParams[:shape]), 1.0)
        end
        Ireco[:,k,j,i] = I
      end
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., acqData.numSlices, acqData.numEchoes, acqData.numCoils)
  return makeAxisArray(Ireco, acqData)
end

"""
  CS Reconstruction using a joint encoding operator for the different echos
  and regularization on the multi-echo data
"""
function reconstruction_multiEcho(acqData::AcquisitionData, recoParams::Dict)

  # encoding operator and trafo into sparse domain
  F = encodingOps2d_multiEcho(acqData,recoParams)

  # sparsifying transform
  recoParams[:multiEcho] = true
  recoParams[:numEchoes] = acqData.numEchoes
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # weight kdata with sampling density
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    weights = vcat(samplingDensity(acqData,recoParams[:shape])...)
  else
    weights = vcat([fill(1.0/sqrt(prod(recoParams[:shape])),size(acqData.kdata[echo,1],2)) for echo=1:acqData.numEchoes]...)
  end
  W = WeightingOp(weights)

  # regularization
  regName = get(recoParams, :regularization, "L1")
  λ = get(recoParams,:λ,0.0)
  normalize = get(recoParams, :normalizeReg, false)

  # reconstruction
  Ireco = zeros(ComplexF64, prod(recoParams[:shape])*acqData.numEchoes, acqData.numCoils, acqData.numSlices)
  solvername = get(recoParams,:solver,"fista")
  for i = 1:acqData.numSlices
    for j = 1:acqData.numCoils
      reg = Regularization(regName, λ; recoParams...)
      if normalize
        RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
      end
      solver = createLinearSolver(solvername, W*F; reg=reg, recoParams...)

      kdata = multiEchoData(acqData, j, i) .* weights  #multiEchoData(acqData2, j, i)
      Ireco[:,j,i] = solve(solver,kdata)
      # TODO circular shutter
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., acqData.numEchoes, acqData.numCoils, acqData.numSlices)
  return makeAxisArray(permutedims(Ireco,[1,2,5,3,4]), acqData)
end

"""
  CS-Sense Reconstruction
"""
function reconstruction_multiCoil(acqData::AcquisitionData, recoParams::Dict)

  # sparsifying transform
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # weight kdata with sampling density
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    weights = samplingDensity(acqData,recoParams[:shape])
  else
    weights = [fill(1.0/sqrt(prod(recoParams[:shape])),size(acqData.kdata[echo,1],2)) for echo=1:acqData.numEchoes]
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  λ = get(recoParams,:λ,0.0)
  normalize = get(recoParams, :normalizeReg, false)

  # solve optimization problem
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), acqData.numSlices, acqData.numEchoes, 1)
  solvername = get(recoParams,:solver,"fista")

  for k = 1:acqData.numSlices
    E = encodingOps2d_parallel(acqData, recoParams; slice=k)
    for j = 1:acqData.numEchoes
      W = WeightingOp(weights[j],acqData.numCoils)
      reg = Regularization(regName, λ; multiEcho=true, recoParams...)
      if normalize
        RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
      end
      solver = createLinearSolver(solvername, W*E[j]; reg=reg, recoParams...)
      kdata = multiCoilData(acqData, j, k) .* repeat(weights[j], acqData.numCoils) #multiCoilData(acqData2, j, k)
      I = solve(solver, kdata)

      if isCircular( trajectory(acqData, j) )
        circularShutter!(reshape(I, recoParams[:shape]), 1.0)
      end
      Ireco[:,j,k] = I
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., acqData.numSlices, acqData.numEchoes, 1)
  return makeAxisArray(Ireco, acqData)
end

function reconstruction_multiCoilMultiEcho(acqData::AcquisitionData, recoParams::Dict)

  # sparsifying transform
  recoParams[:multiEcho] = true
  recoParams[:numEchoes] = acqData.numEchoes
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # rescale data for the symmetrized problem
  densityWeighting = get(recoParams,:densityWeighting,false)
  if densityWeighting
    weights = samplingDensity(acqData,recoParams[:shape])
  else
    weights = vcat([fill(1.0/sqrt(prod(recoParams[:shape])),size(acqData.kdata[echo,1],2)) for echo=1:acqData.numEchoes]...)
  end
  W = WeightingOp(weights,acqData.numCoils)

  # regularization
  regName = get(recoParams, :regularization, "L1")
  λ = get(recoParams,:λ,0.0)
  reg = Regularization(regName, λ; recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
  end

  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), acqData.numEchoes, acqData.numSlices)
  solvername = get(recoParams,:solver,"fista")

  for i = 1:acqData.numSlices
    E = encodingOp_2d_multiEcho_parallel(acqData, recoParams, slice=i)
    solver = createLinearSolver(solvername, W*E; reg=reg, recoParams...)
    kdata = multiCoilMultiEchoData(acqData, i) .* repeat(weights, acqData.numCoils)# multiCoilMultiEchoData(acqData2, i)
    Ireco[:,:,i] = solve(solver, kdata)
  end

  Ireco = reshape( permutedims(Ireco, [1,3,2]), recoParams[:shape]..., acqData.numSlices, acqData.numEchoes )
end
