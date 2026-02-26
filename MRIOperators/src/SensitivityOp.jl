export SensitivityOp

function prod_smap!(y::AbstractVector{T}, smaps::AbstractMatrix{T}, x::AbstractVector{T}, numVox, numChan, numContr=1) where T
  x_ = reshape(x,numVox,numContr)
  y_ = reshape(y,numVox,numChan,numContr)

  @assert size(smaps) == (size(y_,1), size(y_,2))

  @inbounds for i ∈ CartesianIndices(y_)
    y_[i] = x_[i[1],i[3]] * smaps[i[1],i[2]]
  end
  return y
end

function ctprod_smap!(y::AbstractVector{T}, smapsC::AbstractMatrix{T}, x::AbstractVector{T}, numVox, numChan, numContr=1) where T
  x_ = reshape(x,numVox,numChan,numContr)
  y_ = reshape(y,numVox,numContr)

  @assert size(smapsC) == (size(x_,1), size(x_,2))

  y_ .= 0
  @inbounds for i ∈ CartesianIndices(x_)
    y_[i[1],i[3]] += x_[i] * smapsC[i[1],i[2]]
  end
  return y
end


"""
  SensitivityOp(sensMaps::AbstractMatrix{T}, numContr=1)

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`
# Arguments
* `sensMaps`    - sensitivity maps ( 1. dim -> voxels, 2. dim-> coils)
* `numEchoes`   - number of contrasts to which the operator will be applied
"""
struct SensitivityOp{T, O} <: WrappedMRIOperator{T, O}
  op::O
  function SensitivityOp(sensMaps::AbstractMatrix{T}, numContr=1) where T
      numVox, numChan = size(sensMaps)
      sensMapsC = conj.(sensMaps)
      op = LinearOperator{T}(numVox*numContr*numChan, numVox*numContr, false, false,
                           (res,x) -> prod_smap!(res,sensMaps,x,numVox,numChan,numContr),
                           nothing,
                           (res,x) -> ctprod_smap!(res,sensMapsC,x,numVox,numChan,numContr),
                           S = typeof(similar(sensMaps, 0)))
      return new{T, typeof(op)}(op)
  end
end
parent(op::SensitivityOp) = getfield(op, :op)

"""
  SensitivityOp(sensMaps::AbstractArray{T,4}, numContr=1)

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`
# Arguments
* `sensMaps`  - sensitivity maps ( 1.-3. dim -> voxels, 4. dim-> coils)
* `numContr`  - number of contrasts to which the operator will be applied
"""
function SensitivityOp(sensMaps::AbstractArray{T,4}, numContr=1) where T #{T,D}
  sensMaps_mat = reshape(sensMaps, div(length(sensMaps),size(sensMaps,4)),size(sensMaps,4))
  return SensitivityOp(sensMaps_mat, numContr)
end