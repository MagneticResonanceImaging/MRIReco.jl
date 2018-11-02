export reconstruction_direct, reconstruction


function reconstruction_direct(aqData::AcquisitionData, recoParams::Dict)
  numEchoes = aqData.numEchoes
  numCoils = aqData.numCoils
  numSlices = aqData.numSlices
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), numSlices, numEchoes, numCoils)

  aqDataWeighted = weightedData(aqData, recoParams[:shape])

  p = Progress(numSlices*numCoils*numEchoes, 1, "Direct Reconstruction...")

  for i = 1:numSlices
    F = EncodingOp(aqData, recoParams, slice=i)
    for k = 1:numEchoes
      for j = 1:numCoils
        Ireco[:,i,k,j] = adjoint(F[k]) * kData(aqDataWeighted,k,j,i)
        next!(p)
      end
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., numSlices, numEchoes, numCoils)
  return makeAxisArray(Ireco, aqData)
end
