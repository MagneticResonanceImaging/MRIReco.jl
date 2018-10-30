export MapSliceOp

function mapSliceForeward(A, x::Vector, size1::Tuple, dim::Int)
  return vec( mapslices(u->A*u, reshape(x, size1), dims=dim) )
end

function mapSliceBackward(A, y::Vector, size2::Tuple, dim::Int)
  return vec( mapslices(v->adjoint(A)*v, reshape(y, size2), dims=dim) )
end

function MapSliceOp(trafo, dim::Int64, size1::Tuple, size2::Tuple; T=ComplexF64)
  return LinearOperator{T,Function,Nothing,Function}(prod(size2), prod(size1), false, false
            , x->mapSliceForeward(trafo, x, size1, dim)
            , nothing
            , y->mapSliceBackward(trafo, y, size2, dim) )
end
