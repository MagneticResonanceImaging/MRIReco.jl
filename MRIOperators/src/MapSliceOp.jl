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
"""
function MapSliceOp(trafo, dim::Int64, size1::Tuple, size2::Tuple; T=ComplexF64, S = LinearOperators.storage_type(trafo))
  return LinearOperator(prod(size2), prod(size1), false, false
            , (res,x) -> ( res .= mapSliceForeward(trafo, x, size1, dim) )
            , nothing
            , (res,y) -> ( res .= mapSliceBackward(trafo, y, size2, dim) )
            , S = S)
end
