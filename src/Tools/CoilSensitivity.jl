export estimateCoilSensitivities, mergeChannels

"""
  estimateCoilSensitivities(I::AbstractArray{T,5})

Estimates the coil sensitivity based on a reconstruction where the data
from each coil has been reconstructed individually
"""
function estimateCoilSensitivities(I::AbstractArray{T,5}) where T
  numCoils = size(I,5)

  I_sum = sqrt.( sum(abs.(I).^2, dims=5) )

  s = similar(I)
  for i=1:numCoils
    s[:,:,:,:,i] = I[:,:,:,:,i] ./ I_sum
  end

  return s
end

"""
  mergeChannels(I::AbstractArray{T,5})

Merge the channels of a multi-coil reconstruction
"""
mergeChannels(I::AbstractArray{T,5}) where T = sqrt.(sum(abs.(I).^2,dims=5));
