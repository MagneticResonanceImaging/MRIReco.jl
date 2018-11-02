export reconstruction_simple, reconstruction_multiEcho, reconstruction_multiCoil, reconstruction_multiCoilMultiEcho, reconstruction_lowRank

"""
  CS-Sense Reconstruction using sparsity in the wavelet domain
"""
function reconstruction_simple(aqData::AcquisitionData, recoParams::Dict)

  # sparsifying transform
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # weight kdata with sampling density
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    aqData2 = weightedData(aqData, recoParams[:shape])
  else
    aqData2 = aqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, aqData2.kdata)
  end

  # reconstruction
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), aqData.numSlices, aqData.numEchoes, aqData.numCoils)
  solvername = get(recoParams, :solver, "fista")

  for k = 1:aqData.numSlices
    F = EncodingOp(aqData, recoParams, slice=k)
    for j = 1:aqData.numEchoes
      solver = createLinearSolver(solvername, F[j], reg; recoParams...)
      for i = 1:aqData.numCoils
        kdata = kData(aqData2,j,i,k)
        Ireco[:,k,j,i] = solve( solver, kdata )
      end
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., aqData.numSlices, aqData.numEchoes, aqData.numCoils)
  return makeAxisArray(Ireco, aqData)
end

"""
  CS Reconstruction using a joint encoding operator for the different echos
  and regularization on the multi-echo data
"""
function reconstruction_multiEcho(aqData::AcquisitionData, recoParams::Dict)

  # encoding operator and trafo into sparse domain
  F = EncodingOp(aqData,recoParams; multiEcho=true)

  # sparsifying transform
  recoParams[:multiEcho] = true
  recoParams[:numEchoes] = aqData.numEchoes
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # weight kdata with sampling density
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    aqData2 = weightedData(aqData, recoParams[:shape])
  else
    aqData2 = aqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, aqData2.kdata)
  end

  # reconstruction
  Ireco = zeros(ComplexF64, prod(recoParams[:shape])*aqData.numEchoes, aqData.numCoils, aqData.numSlices)
  solvername = get(recoParams,:solver,"fista")
  solver = createLinearSolver(solvername, F, reg; recoParams...)
  for i = 1:aqData.numSlices
    for j = 1:aqData.numCoils
      kdata = multiEchoData(aqData2, j, i)
      Ireco[:,j,i] = solve(solver,kdata)
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., aqData.numEchoes, aqData.numCoils, aqData.numSlices)
  return makeAxisArray(permutedims(Ireco,[1,2,5,3,4]), aqData)
end

"""
  CS-Sense Reconstruction
"""
function reconstruction_multiCoil(aqData::AcquisitionData, recoParams::Dict)

  # sparsifying transform
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # weight kdata with sampling density
  densityWeighting = get(recoParams,:densityWeighting,true)
  if densityWeighting
    aqData2 = weightedData(aqData, recoParams[:shape])
  else
    aqData2 = aqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; multiEcho=true, recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, aqData2.kdata)
  end

  # solve optimization problem
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), aqData.numSlices, aqData.numEchoes, 1)
  solvername = get(recoParams,:solver,"fista")

  for k = 1:aqData.numSlices
    E = EncodingOp(aqData, recoParams; parallel=true, slice=k)
    for j = 1:aqData.numEchoes
      solver = createLinearSolver(solvername, E[j], reg; recoParams...)
      kdata = multiCoilData(aqData2, j, k)
      Ireco[:,j,k] = solve(solver,kdata)
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., aqData.numSlices, aqData.numEchoes, 1)
  return makeAxisArray(Ireco, aqData)
end

function reconstruction_multiCoilMultiEcho(aqData::AcquisitionData, recoParams::Dict)

  # sparsifying transform
  recoParams[:multiEcho] = true
  recoParams[:numEchoes] = aqData.numEchoes
  sparseTrafoName = get(recoParams, :sparseTrafoName, "Wavelet")
  recoParams[:sparseTrafo] = SparseOp(sparseTrafoName; recoParams...)

  # rescale data for the symmetrized problem
  densityWeighting = get(recoParams,:densityWeighting,false)
  if densityWeighting
    aqData2 = weightedData(aqData, recoParams[:shape])
  else
    aqData2 = aqData
  end

  # regularization
  regName = get(recoParams, :regularization, "L1")
  reg = getRegularization(regName; recoParams...)
  normalize = get(recoParams, :normalizeReg, false)
  if normalize
    RegularizedLeastSquares.normalize!(reg, aqData2.kdata)
  end

  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), aqData.numEchoes, aqData.numSlices)
  solvername = get(recoParams,:solver,"fista")

  for i = 1:aqData.numSlices
    E = EncodingOp(aqData, recoParams, parallel=true, multiEcho=true, slice=i)
    solver = createLinearSolver(solvername, E, reg; recoParams...)
    kdata = multiCoilMultiEchoData(aqData2, i)
    Ireco[:,:,i] = solve(solver, kdata)
  end

  Ireco = reshape( permutedims(Ireco, [1,3,2]), recoParams[:shape]..., aqData.numSlices, aqData.numEchoes )
end
