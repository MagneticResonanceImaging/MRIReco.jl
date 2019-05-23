export reconstruction_direct_2d, reconstruction_direct_3d

"""
    reconstruction_direct_2d(acqData::AcquisitionData, reconSize::NTuple{2,Int64}, weights::Vector{Vector{ComplexF64}}, correctionMap::Array{ComplexF64}=ComplexF64[])

Performs a direct Fourier-based image reconstruction of 2d encoded AcquisitionData

input:
  `acqData::AcquisitionData`            - AcquisitionData object
  `reconSize::NTuple{2,Int64}`              - size of image to reconstruct
  `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
  (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects

"""
function reconstruction_direct_2d(acqData::AcquisitionData
                                  , reconSize::NTuple{2,Int64}
                                  , weights::Vector{Vector{ComplexF64}}
                                  , correctionMap::Array{ComplexF64}=ComplexF64[])

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)
  Ireco = zeros(ComplexF64, prod(reconSize), numSl, numContr, numChan)

  p = Progress(numSl*numChan*numContr, 1, "Direct Reconstruction...")

  for i = 1:numSl
    F = encodingOps2d_simple(acqData, reconSize, slice=i, correctionMap=correctionMap)
    for k = 1:numContr
      for j = 1:numChan
        kdata = kData(acqData,k,j,i) .* (weights[k].^2)
        Ireco[:,i,k,j] = adjoint(F[k]) * kdata
        next!(p)
      end
    end
  end

  Ireco = reshape(Ireco, reconSize[1],reconSize[2], numSl, numContr, numChan)
  return makeAxisArray(Ireco, acqData)
end


"""
    reconstruction_direct_3d(acqData::AcquisitionData, reconSize::NTuple{3,Int64}, weights::Vector{Vector{ComplexF64}}, correctionMap::Array{ComplexF64}=ComplexF64[])

Performs a direct Fourier-based image reconstruction of 3d encoded AcquisitionData

input:
  `acqData::AcquisitionData`            - AcquisitionData object
  `reconSize::NTuple{3,Int64}`              - size of image to reconstruct
  `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
  `(correctionMap::Array{ComplexF64})`  - fieldmap for the correction of off-resonance effects

"""
function reconstruction_direct_3d(acqData::AcquisitionData
                                  , reconSize::NTuple{3,Int64}
                                  , weights::Vector{Vector{ComplexF64}}
                                  , correctionMap::Array{ComplexF64}=ComplexF64[])

  numContr, numChan = numContrasts(acqData), numChannels(acqData)
  Ireco = zeros(ComplexF64, prod(reconSize), numContr, numChan)

  p = Progress(numChan*numContr, 1, "Direct Reconstruction...")

  F = encodingOps3d_simple(acqData, reconSize, correctionMap=correctionMap)
  for j = 1:numContr
    for i = 1:numChan
      kdata = kData(acqData,j,i,1) .* (weights[j].^2)
      Ireco[:,j,i] = adjoint(F[j]) * kdata
      next!(p)
    end
  end

  Ireco = reshape(Ireco, reconSize[1], reconSize[2],reconSize[3], numContr, numChan)
  return makeAxisArray(Ireco, acqData)
end

function setupDirectReco(acqData::AcquisitionData, recoParams::Dict)
  reconSize = recoParams[:reconSize]
  weights = samplingDensity(acqData,recoParams[:reconSize])
  # field map
  cmap = get(recoParams, :cmap, ComplexF64[])

  return reconSize, weights, cmap
end
