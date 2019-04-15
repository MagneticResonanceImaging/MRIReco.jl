export reconstruction_direct_2d, reconstruction_direct_3d

function reconstruction_direct_2d(acqData::AcquisitionData, recoParams::Dict)
  numEchoes = acqData.numEchoes
  numCoils = acqData.numCoils
  numSlices = acqData.numSlices
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), numSlices, numEchoes, numCoils)

  # acqDataWeighted = weightedData(acqData, recoParams[:shape])
  weights = samplingDensity(acqData,recoParams[:shape])

  p = Progress(numSlices*numCoils*numEchoes, 1, "Direct Reconstruction...")

  for i = 1:numSlices
    F = encodingOps2d_simple(acqData, recoParams, slice=i)
    for k = 1:numEchoes
      for j = 1:numCoils
        kdata = kData(acqData,k,j,i) .* (weights[k].^2)
        Ireco[:,i,k,j] = adjoint(F[k]) * kdata
        next!(p)
      end
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape][1],recoParams[:shape][2], numSlices, numEchoes, numCoils)
  return makeAxisArray(Ireco, acqData)
end


function reconstruction_direct_3d(acqData::AcquisitionData, recoParams::Dict)
  numEchoes = acqData.numEchoes
  numCoils = acqData.numCoils
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), numEchoes, numCoils)

  # acqDataWeighted = weightedData(acqData, recoParams[:shape])
  weights = samplingDensity(acqData,recoParams[:shape])

  p = Progress(numCoils*numEchoes, 1, "Direct Reconstruction...")

  F = encodingOps3d_simple(acqData, recoParams)
  for j = 1:numEchoes
    for i = 1:numCoils
      kdata = kData(acqData,j,i,1) .* (weights[j].^2)
      Ireco[:,j,i] = adjoint(F[j]) * kdata
      next!(p)
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape][1], recoParams[:shape][2], recoParams[:shape][3], numEchoes, numCoils)
  return makeAxisArray(Ireco, acqData)
end
