export SensitivityOp

#
# normalize sensitivity maps
# ----------------------------------
# 1. dim of sensMaps: voxels
# 2. dim of sensMaps: coils
#
function normalizeSensitivityMap!(sensMaps::Array{T}) where T
  numVoxels, numCoils = size(sensMaps)
  SMap = transpose(sensMaps)
  SMag = zeros(ComplexF64, numVoxels) # calculate the norm over coil-sensitivities for each voxel
  for i = 1:numVoxels
    SMag[i] = norm(SMap[:,i])
  end

  # set values with magnitude zero to infinity, in order to avoid, infinite values
  # when normalizing
  SMag[find(x->x==0.,SMag)] = Inf

  sensMaps[:,:] = sensMaps[:,:]./repeat(SMag,1,numCoils)
end

function sensOpForeward(sensMaps::Matrix, x::Vector{T}, numEchoes::Int=1) where T
  return  vec( repeat(sensMaps,numEchoes,1).*x )
end

function sensOpBackward(sensMaps::Matrix, y::Vector{T}, numEchoes::Int=1) where T
  y = reshape( y, size(sensMaps,1)*numEchoes, size(sensMaps,2) )
  return vec( sum( repeat(conj(sensMaps),numEchoes,1).*y, dims=2 ) )
end

function SensitivityOp(sensMaps::Matrix{T}, numEchoes::Int=1) where T
  return LinearOperator{T}(length(sensMaps)*numEchoes, size(sensMaps,1)*numEchoes, false, false
      , x->sensOpForeward(sensMaps, x, numEchoes)
      , nothing
      , y-> sensOpBackward(sensMaps, y, numEchoes) )
end

function SensitivityOp(sensMaps::Array{T,D}, numEchoes::Int=1) where {T,D}
  return SensitivityOp(reshape(sensMaps, div(length(sensMaps),size(sensMaps,D)),
                                         size(sensMaps,D)), numEchoes)
end
