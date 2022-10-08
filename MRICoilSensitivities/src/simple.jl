
"""
    `s = estimateCoilSensitivities(I::AbstractArray{T,6})`

Estimates the coil sensitivity based on a reconstruction where the data
from each coil has been reconstructed individually.
Returns a 5D array.
"""
function estimateCoilSensitivities(I::AbstractArray{T,6}, thresh = 1.e-2) where {T}
  nx, ny, nz, ne, numChan = size(I)

  I_sum = sqrt.(sum(abs.(I) .^ 2, dims = 5)) .+ eps()
  I_max = maximum(abs.(I_sum))
  msk = zeros(size(I_sum))
  msk[findall(x -> x > thresh * I_max, I_sum)] .= 1

  s = zeros(eltype(I), size(I))
  for i = 1:numChan
    s[:, :, :, :, i] = msk .* I[:, :, :, :, i] ./ I_sum
  end

  return s
end


"""
    `I4 = mergeChannels(I::AbstractArray{T,6})`

Merge the channels of a multi-coil reconstruction.
Returns a 6D array.
"""
mergeChannels(I::AbstractArray{T,6}) where {T} = sqrt.(sum(abs.(I) .^ 2, dims = 5))
