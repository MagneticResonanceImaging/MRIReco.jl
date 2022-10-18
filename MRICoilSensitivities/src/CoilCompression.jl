"""
return SVD-based coil compression matrix for `numVC` virtual coils
"""
function geometricCCMat(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2)) where T <: Complex
  return svd(kdata).Vt[:, 1:numVC]
end

"""
perform SVD-based coil compression for the given `kdata` and `smaps`
"""
function geometricCoilCompression(kdata::Matrix{T}, smaps::Array{T,4}, numVC::Int64 = size(kdata, 2)) where T <: Complex
  nx, ny, nz, nc = size(smaps)
  usv = svd(kdata)
  ccMat = usv.Vt[:, 1:numVC]
  kdataCC = kdata * ccMat
  smapsCC = zeros(T, nx, ny, nz, numVC)
  for j = 1:nz, i = 1:ny
    smapsCC[:, i, j, :] .= smaps[:, i, j, :] * ccMat
  end

  return kdataCC, smapsCC
end

"""
perform SVD-based coil compression for the given `acqData` and `smaps`
"""
function geometricCoilCompression(acqData::AcquisitionData{T,2}, smaps::Array{Complex{T},4}, numVC::Int64 = size(smaps, 4)) where T
  nx, ny, nz, nc = size(smaps)
  acqDataCC = deepcopy(acqData)
  smapsCC = zeros(Complex{T}, nx, ny, nz, numVC)
  for sl = 1:numSlices(acqData)
    # use first echo and first repetition to determine compression matrix
    usv = svd(acqData.kdata[1, sl, 1])
    ccMat = usv.Vt[:, 1:numVC]
    # compress slice of the smaps
    for i = 1:ny
      smapsCC[:, i, sl, :] .= smaps[:, i, sl, :] * ccMat
    end
    # compress kdata-slice
    for rep = 1:numRepetitions(acqData), contr = 1:numContrasts(acqData)
      acqDataCC.kdata[contr, sl, rep] = acqData.kdata[contr, sl, rep] * ccMat
    end
  end
  return acqDataCC, smapsCC
end

function geometricCoilCompression(acqData::AcquisitionData{T,3}, smaps::Array{Complex{T},4}, numVC::Int64 = size(smaps, 4)) where T
  nc = size(smaps, 4)
  N = size(smaps)[1:3]
  acqDataCC = deepcopy(acqData)
  smapsCC = zeros(Complex{T}, N..., numVC)
  
  # use first echo and first repetition to determine compression matrix
  usv = svd(acqData.kdata[1, 1, 1])
  ccMat = usv.Vt[:, 1:numVC]

  # compress slice of the smaps
  for n in CartesianIndices(N)
    smapsCC[n, :] .= transpose(ccMat) * smaps[n, :]
  end

  # compress kdata-slice
  for rep = 1:numRepetitions(acqData), contr = 1:numContrasts(acqData)
    acqDataCC.kdata[contr, 1, rep] = acqData.kdata[contr, 1, rep] * ccMat
  end
  
  return acqDataCC, smapsCC
end
