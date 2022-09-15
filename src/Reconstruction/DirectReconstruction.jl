export reconstruction_direct_2d, reconstruction_direct_3d

"""
    reconstruction_direct(acqData::AcquisitionData, reconSize::NTuple{D,Int64}, weights::Vector{Vector{Complex{<:AbstractFloat}}}, correctionMap::Array{Complex{<:AbstractFloat}}=Complex{<:AbstractFloat}[])

Performs a direct Fourier-based image reconstruction of AcquisitionData

input:
  `acqData::AcquisitionData`            - AcquisitionData object
  `reconSize::NTuple{D,Int64}`              - size of image to reconstruct
  `weights::Vector{Vector{Complex{<:AbstractFloat}}}` - sampling density of the trajectories in acqData
  (`correctionMap::Array{Complex{<:AbstractFloat}}`)  - fieldmap for the correction of off-resonance effects

"""
function reconstruction_direct(acqData::AcquisitionData{T}
                                  , reconSize::NTuple{D,Int64}
                                  , weights::Vector{Vector{Complex{T}}}
                                  , correctionMap::Array{Complex{T}}=Complex{T}[]) where {D,T}

  encDims = ndims(trajectory(acqData))
  if encDims!=D
    error("reco-dimensionality $D and encoding-dimensionality $(encDims) do not match")
  end

  numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)
  Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numChan, numRep)

  p = Progress(numSl*numChan*numContr*numRep, 1, "Direct Reconstruction...")

  for i = 1:numSl
    F = encodingOps_simple(acqData, reconSize, slice=i, correctionMap=correctionMap)
    for k = 1:numContr
      for j = 1:numChan
        for l = 1:numRep
            kdata = kData(acqData,k,j,i,rep=l) .* (weights[k].^2)
            I = adjoint(F[k]) * kdata

            if isCircular( trajectory(acqData, k) )
            circularShutter!(reshape(I, reconSize), 1.0)
            end
            Ireco[:,i,k,j,l] = I

            next!(p)
        end
      end
    end
  end

  Ireco = reshape(Ireco, volumeSize(reconSize, numSl)..., numContr, numChan,numRep)
  # if D==2
  #   # 2d reconstruction
  #   Ireco = reshape(Ireco, reconSize[1], reconSize[2], numSl, numContr, numChan)
  # else
  #   # 3d reconstruction
  #   Ireco = reshape(Ireco, reconSize[1], reconSize[2], reconSize[3], numContr, numChan)
  # end

  return makeAxisArray(Ireco, acqData)
end
