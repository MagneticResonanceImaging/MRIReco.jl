export MapSliceOp

function mapSliceForeward(A, x::AbstractVector, size1::Tuple, dim::Int)
  return vec( mapslices(u->A*u, reshape(x, size1), dims=dim) )
end

function mapSliceBackward(A, y::AbstractVector, size2::Tuple, dim::Int)
  return vec( mapslices(v->adjoint(A)*v, reshape(y, size2), dims=dim) )
end

"""
    MapSliceOp(trafo, dim::Int64, size1::Tuple, size2::Tuple; T=ComplexF64)

generates an operator that applies `trafo` to each slice of dimension `dim` of an array

# Arguments
* `trafo`           - transformation to apply
* `size1`           - size of the Array to be transformed
* `size2`           - size of the resulting Array of applying trafo to all slices
* (`T=ComplexF64`)  - type of the transformation

# Moving to another array backend

    Adapt.adapt(CuArray, op::MapSliceOp)

returns a new `MapSliceOp` with a new parent operator, built from `Adapt.adapt(CuArray, trafo)` and the same
`dim`, `size1`, `size2` and `T`. An operator `trafo` is adapted by its own Adapt rule, such as those of
LinearOperatorCollection; an operator without one is kept unchanged. On the same backend, an array `trafo` may be
shared with `op` and is only read; an operator `trafo` kept unchanged is shared with `op`, including any internal
state it updates when applied. `op` itself is not modified.
`SensitivityOp`, `SubspaceOp`, `MapSliceOp` and `ExplicitOp` are the MRI operators that support `Adapt.adapt` so far.
"""
struct MapSliceOp{T, O, A, S1 <: Tuple, S2 <: Tuple} <: WrappedMRIOperator{T, O}
  op::O
  trafo::A
  dim::Int
  size1::S1
  size2::S2
  function MapSliceOp(trafo, dim::Int64, size1::Tuple, size2::Tuple; T=ComplexF64, S = LinearOperators.storage_type(trafo))
    op = LinearOperator{T}(prod(size2), prod(size1), false, false
              , (res,x) -> ( res .= mapSliceForeward(trafo, x, size1, dim) )
              , nothing
              , (res,y) -> ( res .= mapSliceBackward(trafo, y, size2, dim) )
              , S = S)
    return new{T, typeof(op), typeof(trafo), typeof(size1), typeof(size2)}(op, trafo, dim, size1, size2)
  end
end
parent(op::MapSliceOp) = getfield(op, :op)

Adapt.adapt_structure(to, op::MapSliceOp{T}) where T =
  MapSliceOp(adapt(to, getfield(op, :trafo)), getfield(op, :dim), getfield(op, :size1), getfield(op, :size2); T,
             S = typeof(adapt(to, LinearOperators.storage_type(op)(undef, 0))))
