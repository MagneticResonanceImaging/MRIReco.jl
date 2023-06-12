export applySubspace
"""
"""
function applySubspace(I::AbstractArray{Complex{T},6}, basis::AbstractMatrix{Complex{T}}) where {T}
  sx,sy,sz,numEchoes,numCoils,numReps = size(I)
  I = reshape(I,sx*sy*sz,numEchoes,numCoils,numReps)

  I2 = zeros(Complex{T},sx*sy*sz,size(basis,1),numCoils,numReps)
  for i in 1:numCoils, j in 1:numReps
    I2[:,:,i,j] .= I[:,:,i,j] * basis'
  end

  I2 = reshape(I2,sx,sy,sz,:,numCoils,numReps)
  return I2
end

function applySubspace(I::AxisArray{Complex{T},6}, basis::AbstractMatrix{Complex{T}}) where {T}
  axesVec = AxisArrays.axes(I)
  I2 = applySubspace(I.data,basis)

  I2 = AxisArray(I2,
  axesVec[1:3]...,
  Axis{:echos}(1:size(I2,4)),
  axesVec[5:6]...)

  return I2
end