"""
return SVD-based coil compression matrix for `numVC` virtual coils
"""
function geometricCCMat(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2)) where T <: Complex
  return svd(kdata).Vt[:, 1:numVC]
end

"""
perform SVD-based coil compression for the given `kdata` and `smaps`
"""
function geometricCC(kdata::Matrix{T}, smaps::Array{T,4}, numVC::Int64 = size(kdata, 2)) where T <: Complex
  nx, ny, nz, nc = size(smaps)
  usv = svd(kdata)
  kdataCC = kdata * usv.Vt[:, 1:numVC]
  smapsCC = zeros(T, nx, ny, nz, numVC)
  for j = 1:nz, i = 1:ny
    smapsCC[:, i, j, :] .= smaps[:, i, j, :] * usv.Vt[:, 1:numVC]
  end

  return kdataCC, smapsCC
end

"""
perform SVD-based coil compression for the given 2d-encoded `acqData` and `smaps`
"""
function geometricCC_2d(acqData::AcquisitionData{T}, smaps::Array{Complex{T},4}, numVC::Int64 = size(smaps, 4)) where T
  nx, ny, nz, nc = size(smaps)
  acqDataCC = deepcopy(acqData)
  smapsCC = zeros(Complex{T}, nx, ny, nz, numVC)
  for sl = 1:numSlices(acqData)
    # use first echo and first repetition to determine compression matrix
    ccMat = geometricCCMat(acqData.kdata[1, sl, 1], numVC)
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
