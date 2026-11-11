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

# Moving to another array backend

    Adapt.adapt(CuArray, op::SubspaceOp)

returns a new `SubspaceOp` with a new parent operator, built from `Adapt.adapt(CuArray, basis)`.
Its storage type `S` is derived from the adapted basis. On the same backend, the basis may be shared with
`op` and is only read. `op` itself is not modified.
`SensitivityOp`, `SubspaceOp` and `FieldmapNFFTOp` are the MRI operators that support `Adapt.adapt` so far.
"""
struct SubspaceOp{T, O, M <: AbstractMatrix{T}, D} <: WrappedMRIOperator{T, O}
  op::O
  basis::M
  shape::NTuple{D, Int}
  numContr::Int
  function SubspaceOp(basis::AbstractMatrix{T},shape::NTuple{D,Int64},numContr, S = LinearOperators.storage_type(basis)) where {T,D}
      numVox = prod(shape)
      numBasis = size(basis,2)
      op = LinearOperator{T}(numVox*numContr, numVox*numBasis, false, false,
                           (res,x) -> prod_subspace!(res,basis,x,numVox,numContr,numBasis),
                           nothing,
                           (res,x) -> ctprod_subspace!(res,basis,x,numVox,numContr,numBasis),
                           S = S)
      return new{T, typeof(op), typeof(basis), D}(op, basis, shape, numContr)
  end
end
parent(op::SubspaceOp) = getfield(op, :op)

Adapt.adapt_structure(to, op::SubspaceOp) =
  SubspaceOp(adapt(to, getfield(op, :basis)), getfield(op, :shape), getfield(op, :numContr))