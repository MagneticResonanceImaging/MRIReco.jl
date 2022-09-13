"""
  SensitivityOp(sensMaps::Matrix{ComplexF64}, numEchoes::Int=1)

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`

# Arguments
* `sensMaps::Matrix{ComplexF64}`  - sensitivity maps ( 1. dim -> voxels, 2. dim-> coils)
* `numEchoes`                     - number of contrasts to which the opetaor will be applied
"""
function SensitivityOp(sensMaps::Matrix{T}, numContr::Int=1) where T
  sensOps = [opDiagonal(sensMaps[:,c]) for contr=1:numContr, c=1:size(sensMaps,2)]
  return vcat( sensOps...)
end

"""
  SensitivityOp(sensMaps::Array{T,4}, numContr::Int=1) where T

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`

# Arguments
* `sensMaps::Array{T,4}`  - sensitivity maps ( 1.-3. dim -> voxels, 4. dim-> coils)
* `numContr`             - number of contrasts to which the opetaor will be applied
"""
function SensitivityOp(sensMaps::Array{T,4}, numContr::Int=1) where T #{T,D}
sensMaps_mat = reshape(sensMaps, div(length(sensMaps),size(sensMaps,4)),size(sensMaps,4))
return SensitivityOp(sensMaps_mat, numContr)
end