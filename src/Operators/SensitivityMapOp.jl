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

"""
    SensitivityOp(sensMaps::Matrix{ComplexF64}, numEchoes::Int=1)

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`

# Arguments
* `sensMaps::Matrix{ComplexF64}`  - sensitivity maps ( 1. dim -> voxels, 2. dim-> coils)
* `numEchoes`                     - number of contrasts to which the opetaor will be applied
"""
function SensitivityOp(sensMaps::Matrix{ComplexF64}, numEchoes::Int=1)
  sensMaps_cat = repeat(sensMaps,numEchoes,1)
  return vcat( [opDiagonal(sensMaps_cat[:,coil]) for coil=1:size(sensMaps,2)]... )
end

"""
    SensitivityOp(sensMaps::Array{T,4}, numEchoes::Int=1) where T

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`

# Arguments
* `sensMaps::Array{T,4}`  - sensitivity maps ( 1.-3. dim -> voxels, 4. dim-> coils)
* `numEchoes`             - number of contrasts to which the opetaor will be applied
"""
function SensitivityOp(sensMaps::Array{T,4}, numEchoes::Int=1) where T #{T,D}
  sensMaps_mat = reshape(sensMaps, div(length(sensMaps),size(sensMaps,4)),size(sensMaps,4))
  return SensitivityOp(sensMaps_mat, numEchoes)
end
