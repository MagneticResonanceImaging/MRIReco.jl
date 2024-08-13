export SubspaceOp

function prod_subspace!(y::AbstractVector{T}, basis::AbstractMatrix{T}, x::AbstractVector{T}, numVox, numContr, numBasis) where T
  x_ = reshape(x,numVox,numBasis)
  y_ = reshape(y,numVox,numContr)
  
  y_ .= x_ * basis'
  return y
end

function ctprod_subspace!(y::AbstractVector{T},basis::AbstractMatrix{T}, x::AbstractVector{T}, numVox, numContr, numBasis) where T
  x_ = reshape(x,numVox,numContr)
  y_ = reshape(y,numVox,numBasis)

  y_ .= x_ * basis
  return y
end


"""
  SubspaceOp(basis::AbstractMatrix{T},shape::NTuple{D,Int64},numContr) where {T,D}

basis correspond to svd_object.V cropped to a certain level

basis type should correspond to the type of the rawdata
"""
function SubspaceOp(basis::AbstractMatrix{T},shape::NTuple{D,Int64},numContr, S = LinearOperators.storage_type(basis)) where {T,D}
    numVox = prod(shape)
    numBasis = size(basis,2)

    return LinearOperator{T}(numVox*numContr, numVox*numBasis, false, false,
                         (res,x) -> prod_subspace!(res,basis,x,numVox,numContr,numBasis),
                         nothing,
                         (res,x) -> ctprod_subspace!(res,basis,x,numVox,numContr,numBasis),
                         S = S)
end


