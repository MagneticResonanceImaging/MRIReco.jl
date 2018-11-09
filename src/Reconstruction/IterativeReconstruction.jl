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
    acqData2 = weightedData(acqData, recoParams[:shape])
  else
    acqData2 = acqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
  end

  # reconstruction
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), acqData.numSlices, acqData.numEchoes, acqData.numCoils)
  solvername = get(recoParams, :solver, "fista")

  for k = 1:acqData.numSlices
    F = EncodingOp2d(acqData, recoParams, slice=k)
    for j = 1:acqData.numEchoes
      solver = createLinearSolver(solvername, F[j], reg; recoParams...)
      for i = 1:acqData.numCoils
        kdata = kData(acqData2,j,i,k)
        Ireco[:,k,j,i] = solve( solver, kdata )
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
  F = EncodingOp2d(acqData,recoParams; multiEcho=true)

  # sparsifying transform
  recoParams[:multiEcho] = true
  recoParams[:numEchoes] = acqData.numEchoes
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # weight kdata with sampling density
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    acqData2 = weightedData(acqData, recoParams[:shape])
  else
    acqData2 = acqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
  end

  # reconstruction
  Ireco = zeros(ComplexF64, prod(recoParams[:shape])*acqData.numEchoes, acqData.numCoils, acqData.numSlices)
  solvername = get(recoParams,:solver,"fista")
  solver = createLinearSolver(solvername, F, reg; recoParams...)
  for i = 1:acqData.numSlices
    for j = 1:acqData.numCoils
      kdata = multiEchoData(acqData2, j, i)
      Ireco[:,j,i] = solve(solver,kdata)
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
    acqData2 = weightedData(acqData, recoParams[:shape])
  else
    acqData2 = acqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; multiEcho=true, recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
  end

  # solve optimization problem
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), acqData.numSlices, acqData.numEchoes, 1)
  solvername = get(recoParams,:solver,"fista")

  for k = 1:acqData.numSlices
    E = EncodingOp2d(acqData, recoParams; parallel=true, slice=k)
    for j = 1:acqData.numEchoes
      solver = createLinearSolver(solvername, E[j], reg; recoParams...)
      kdata = multiCoilData(acqData2, j, k)
      Ireco[:,j,k] = solve(solver,kdata)
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
    acqData2 = weightedData(acqData, recoParams[:shape])
  else
    acqData2 = acqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, acqData2.kdata)
  end

  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), acqData.numEchoes, acqData.numSlices)
  solvername = get(recoParams,:solver,"fista")

  for i = 1:acqData.numSlices
    E = EncodingOp2d(acqData, recoParams, parallel=true, multiEcho=true, slice=i)
    solver = createLinearSolver(solvername, E, reg; recoParams...)
    kdata = multiCoilMultiEchoData(acqData2, i)
    Ireco[:,:,i] = solve(solver, kdata)
  end

  Ireco = reshape( permutedims(Ireco, [1,3,2]), recoParams[:shape]..., acqData.numSlices, acqData.numEchoes )
end
