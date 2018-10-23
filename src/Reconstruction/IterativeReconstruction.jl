export reconstruction_simple, reconstruction_multiEcho, reconstruction_multiCoil, reconstruction_multiCoilMultiEcho, reconstruction_lowRank

"""
  CS-Sense Reconstruction using sparsity in the wavelet domain
"""
function reconstruction_simple(aqData::AquisitionData, recoParams::Dict)

  # operators
  F = EncodingOp(aqData, recoParams;numEchoes=aqData.numEchoes)

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
    LinearSolver.normalize!(reg, aqData2.kdata)
  end

  # reconstruction
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), aqData.numSlices, aqData.numEchoes, aqData.numCoils)
  solvername = get(recoParams,:solver,"fista")
  solver = [createLinearSolver(solvername, F[i], reg; recoParams...) for i=1:aqData.numEchoes]
  for i = 1:aqData.numCoils
    for j = 1:aqData.numEchoes
      for k = 1:aqData.numSlices
        kdata = kData(aqData2,k,j,i)
        Ireco[:,k,j,i] = solve( solver[j], kdata )
      end
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., aqData.numSlices, aqData.numEchoes, aqData.numCoils)
end

"""
  CS Reconstruction using a joint encoding operator for the different echos
  and regularization on the multi-echo data
"""
function reconstruction_multiEcho(aqData::AquisitionData, recoParams::Dict)

  # encoding operator and trafo into sparse domain
  F = EncodingOp(aqData,recoParams; numEchoes=aqData.numEchoes, multiEcho=true)

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
    LinearSolver.normalize!(reg, aqData2.kdata)
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
  return permutedims(Ireco,[1,2,5,3,4])
end

"""
  CS-Sense Reconstruction
"""
function reconstruction_multiCoil(aqData::AquisitionData, recoParams::Dict)

  # encoding operator and trafo into sparse domain
  E = EncodingOp(aqData, recoParams; numEchoes=aqData.numEchoes, parallel=true)

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
    LinearSolver.normalize!(reg, aqData2.kdata)
  end

  # solve optimization problem
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), aqData.numSlices, aqData.numEchoes, 1)
  solvername = get(recoParams,:solver,"fista")
  solver = [ createLinearSolver(solvername, E[j], reg; recoParams...) for j=1:aqData.numEchoes]
  for i = 1:aqData.numEchoes
    for j = 1:aqData.numSlices
      kdata = multiCoilData(aqData2, j, i)
      Ireco[:,j,i] = solve(solver[i],kdata)
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., aqData.numSlices, aqData.numEchoes, 1)
end

function reconstruction_multiCoilMultiEcho(aqData::AquisitionData, recoParams::Dict)

  # signal encoding transformation
  E = EncodingOp(aqData, recoParams, numEchoes=aqData.numEchoes, parallel=true, multiEcho=true)

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
    LinearSolver.normalize!(reg, aqData2.kdata)
  end

  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), aqData.numEchoes, aqData.numSlices)
  solvername = get(recoParams,:solver,"fista")
  solver = createLinearSolver(solvername, E, reg; recoParams...)
  for i = 1:aqData.numSlices
    kdata = multiCoilMultiEchoData(aqData2,i)
    Ireco[:,:,i] = solve(solver, kdata)
  end

  Ireco = reshape( permutedims(Ireco, [1,3,2]), recoParams[:shape]..., aqData.numSlices, aqData.numEchoes )
end
