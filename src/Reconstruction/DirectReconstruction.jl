export reconstruction_direct_2d, reconstruction_direct_3d

"""
    reconstruction_direct(acqData::AcquisitionData, reconSize::NTuple{D,Int64}, weights::Vector{Vector{ComplexF64}}, correctionMap::Array{ComplexF64}=ComplexF64[])

Performs a direct Fourier-based image reconstruction of AcquisitionData

input:
  `acqData::AcquisitionData`            - AcquisitionData object
  `reconSize::NTuple{D,Int64}`              - size of image to reconstruct
  `weights::Vector{Vector{ComplexF64}}` - sampling density of the trajectories in acqData
  (`correctionMap::Array{ComplexF64}`)  - fieldmap for the correction of off-resonance effects

"""
function reconstruction_direct(acqData::AcquisitionData
                                  , reconSize::NTuple{D,Int64}
                                  , weights::Vector{Vector{ComplexF64}}
                                  , correctionMap::Array{ComplexF64}=ComplexF64[]) where D

  encDims = dims(trajectory(acqData))
  if encDims!=D
    error("reco-dimensionality $D and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)
  Ireco = zeros(ComplexF64, prod(reconSize), numSl, numContr, numChan)

  p = Progress(numSl*numChan*numContr, 1, "Direct Reconstruction...")

  for i = 1:numSl
    F = encodingOps_simple(acqData, reconSize, slice=i, correctionMap=correctionMap)
    for k = 1:numContr
      for j = 1:numChan
        kdata = kData(acqData,k,j,i) .* (weights[k].^2)
        I = adjoint(F[k]) * kdata

        if isCircular( trajectory(acqData, k) )
          circularShutter!(reshape(I, reconSize), 1.0)
        end
        Ireco[:,i,k,j] = I

        next!(p)
      end
    end
  end

  if D==2
    # 2d reconstruction
    Ireco = reshape(Ireco, reconSize[1], reconSize[2], numSl, numContr, numChan)
  else
    # 3d reconstruction
    Ireco = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numContr, numChan)
  end

  return makeAxisArray(Ireco, acqData)
end
