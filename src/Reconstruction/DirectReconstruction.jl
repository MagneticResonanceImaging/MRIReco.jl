export reconstruction_direct_2d, reconstruction_direct_3d

function reconstruction_direct_2d(acqData::AcquisitionData
                                  , shape::NTuple{2,Int64}
                                  , weights::Vector{Vector{ComplexF64}}
                                  , correctionMap::Array{ComplexF64}=ComplexF64[])
  numEchoes = acqData.numEchoes
  numCoils = acqData.numCoils
  numSlices = acqData.numSlices
  Ireco = zeros(ComplexF64, prod(shape), numSlices, numEchoes, numCoils)

  p = Progress(numSlices*numCoils*numEchoes, 1, "Direct Reconstruction...")

  for i = 1:numSlices
    F = encodingOps2d_simple(acqData, shape, slice=i, correctionMap=correctionMap)
    for k = 1:numEchoes
      for j = 1:numCoils
        kdata = kData(acqData,k,j,i) .* (weights[k].^2)
        Ireco[:,i,k,j] = adjoint(F[k]) * kdata
        next!(p)
      end
    end
  end

  Ireco = reshape(Ireco, shape[1],shape[2], numSlices, numEchoes, numCoils)
  return makeAxisArray(Ireco, acqData)
end


# function reconstruction_direct_3d(acqData::AcquisitionData, recoParams::Dict)
function reconstruction_direct_3d(acqData::AcquisitionData
                                  , shape::NTuple{3,Int64}
                                  , weights::Vector{Vector{ComplexF64}}
                                  , correctionMap::Array{ComplexF64}=ComplexF64[])
  numEchoes = acqData.numEchoes
  numCoils = acqData.numCoils
  Ireco = zeros(ComplexF64, prod(shape), numEchoes, numCoils)

  p = Progress(numCoils*numEchoes, 1, "Direct Reconstruction...")

  F = encodingOps3d_simple(acqData, shape, correctionMap=correctionMap)
  for j = 1:numEchoes
    for i = 1:numCoils
      kdata = kData(acqData,j,i,1) .* (weights[j].^2)
      Ireco[:,j,i] = adjoint(F[j]) * kdata
      next!(p)
    end
  end

  Ireco = reshape(Ireco, shape[1], shape[2],shape[3], numEchoes, numCoils)
  return makeAxisArray(Ireco, acqData)
end

function setupDirectReco(acqData::AcquisitionData, recoParams::Dict)
  shape = recoParams[:shape]
  weights = samplingDensity(acqData,recoParams[:shape])
  # field map
  cmap = get(recoParams, :cmap, ComplexF64[])

  return shape, weights, cmap
end
