export reconstruction_direct, reconstruction


function reconstruction_direct(aqData::AquisitionData, recoParams::Dict)
  numEchoes = aqData.numEchoes
  numCoils = aqData.numCoils
  numSlices = aqData.numSlices
  Ireco = zeros(ComplexF64, prod(recoParams[:shape]), numEchoes, numCoils, numSlices)

  F = EncodingOp(aqData, recoParams; numEchoes=numEchoes)

  aqDataWeighted = weightedData(aqData, recoParams[:shape])

  p = Progress(numSlices*numCoils*numEchoes, 1, "Direct Reconstruction...")
  for i = 1:numSlices
    for j= 1:numCoils
      for k = 1:numEchoes
        Ireco[:,k,j,i] =   adjoint(F[k]) * kData(aqDataWeighted,k,j,i)
        next!(p)
      end
    end
  end

  Ireco = reshape(Ireco, recoParams[:shape]..., numEchoes, numCoils, numSlices)
end
