export reconstruction_direct, reconstruction


function reconstruction_direct(acqData::AcquisitionData, recoParams::Dict)
  numEchoes = acqData.numEchoes
  numCoils = acqData.numCoils
  numSlices = acqData.numSlices
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), numSlices, numEchoes, numCoils)

  acqDataWeighted = weightedData(acqData, recoParams[:shape])

  p = Progress(numSlices*numCoils*numEchoes, 1, "Direct Reconstruction...")

  for i = 1:numSlices
    F = EncodingOp(acqData, recoParams, slice=i)
    for k = 1:numEchoes
      for j = 1:numCoils
        Ireco[:,i,k,j] = adjoint(F[k]) * kData(acqDataWeighted,k,j,i)
        next!(p)
      end
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., numSlices, numEchoes, numCoils)
  return makeAxisArray(Ireco, acqData)
end
