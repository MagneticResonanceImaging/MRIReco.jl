export SensitivityOp

#
# normalize sensitivity maps
# ----------------------------------
# 1. dim of sensMaps: voxels
# 2. dim of sensMaps: coils
#
# FIXME: assumes only one slice
function normalizeSensitivityMap!(sensMaps::Matrix{ComplexF64})
  numVoxels, numCoils = size(sensMaps)
  SMap = transpose(reshape(sensMaps,numVoxels,numCoils))
  SMag = zeros(ComplexF64, numVoxels) # calculate the norm over coil-sensitivities for each voxel
  for i = 1:numVoxels
    SMag[i] = norm(SMap[:,i])
  end

  # set values with magnitude zero to infinity, in order to avoid, infinite values
  # when normalizing
  SMag[find(x->x==0.,SMag)] = Inf

  sensMaps[:,:] = sensMaps[:,:]./repeat(SMag,1,numCoils)
end

function sensOpForeward(sensMaps::Matrix{ComplexF64}, x::Vector{T}, numEchoes::Int=1) where T
  return  vec( repeat(sensMaps,numEchoes,1).*x )
end

function sensOpBackward(sensMaps::Matrix{ComplexF64}, y::Vector{ComplexF64}, numEchoes::Int=1)
  y = reshape( y, size(sensMaps,1)*numEchoes, size(sensMaps,2) )
  return vec( sum( repeat(conj(sensMaps),numEchoes,1).*y, dims=2 ) )
end

function SensitivityOp(sensMaps::Matrix{ComplexF64}, numEchoes::Int=1)
  return LinearOperator{ComplexF64,Function,Nothing,Function}(length(sensMaps)*numEchoes, size(sensMaps,1)*numEchoes, false, false
      , x->sensOpForeward(sensMaps, x, numEchoes)
      , nothing
      , y-> sensOpBackward(sensMaps, y, numEchoes) )
end

function SensitivityOp(sensMaps::Array{T,4}, numEchoes::Int=1) where T #{T,D}
  return SensitivityOp(reshape(sensMaps, div(length(sensMaps),size(sensMaps,4)),
                                         size(sensMaps,4)), numEchoes)
end
