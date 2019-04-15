import Base: hcat, vcat, \
export hcat, vcat, \, diagOp, A_mul_B

include("NFFTOp.jl")
include("ExplicitOp.jl")
include("SensitivityMapOp.jl")
include("SamplingOp.jl")
include("MapSliceOp.jl")
include("FieldmapNFFTOp.jl")
include("EncodingOp.jl")
include("SparseOp.jl")

# FIXME: Weighting Op should be removed, as it is contained on the
# master branch of RegularizedLeastSquares.jl
include("WeightingOp.jl")

#
# horizontally concatenate linear operators n times
#
function hcat(A::AbstractLinearOperator, n::Int)
  op = A
  for i = 2:n
    op = LinearOperators.hcat(op, A)
  end
  return op
end

#
# vertically concatenate linear operators n times
#
function vcat(A::AbstractLinearOperator, n::Int)
  op = A
  for i = 2:n
    op = LinearOperators.vcat(op, A)
  end
  return op
end

function diagOpProd(x::Vector{T}, nrow::Int, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, nrow)
  xIdx=1
  yIdx=1
  for i=1:length(ops)
    y[yIdx:yIdx+ops[i].nrow-1] = ops[i]*x[xIdx:xIdx+ops[i].ncol-1]
    xIdx += ops[i].ncol
    yIdx += ops[i].nrow
  end
  return y
end

function diagOpTProd(x::Vector{T}, ncol::Int, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, ncol)
  xIdx=1
  yIdx=1
  for i=1:length(ops)
    # y[yIdx:yIdx+ops[i].ncol-1] = (ops[i].')*x[xIdx:xIdx+ops[i].nrow-1]
    y[yIdx:yIdx+ops[i].ncol-1] = transpose(ops[i])*x[xIdx:xIdx+ops[i].nrow-1]
    xIdx += ops[i].nrow
    yIdx += ops[i].ncol
  end
  return y
end

function diagOpCTProd(x::Vector{T}, ncol::Int, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, ncol)
  xIdx=1
  yIdx=1
  for i=1:length(ops)
    y[yIdx:yIdx+ops[i].ncol-1] = (ops[i])'*x[xIdx:xIdx+ops[i].nrow-1]
    xIdx += ops[i].nrow
    yIdx += ops[i].ncol
  end
  return y
end

#
# create a bloc-diagonal operator out of smaller operators
#
function diagOp(ops :: AbstractLinearOperator...)
  nrow=0
  ncol=0
  S = eltype(ops[1])
  for i = 1:length(ops)
    nrow += ops[i].nrow
    ncol += ops[i].ncol
    S = promote_type(S, eltype(ops[i]))
  end

  return LinearOperator{S,Function,Function,Function}( nrow, ncol, false, false
                            , x->diagOpProd(x,nrow,ops...)
                            , y->diagOpTProd(y,ncol,ops...)
                            , y->diagOpCTProd(y,ncol,ops...) )
end

# implement A_mul_B for the product
A_mul_B(A::AbstractLinearOperator, x::Vector) = A*x

#
# use hermitian conjugate as an estimate for the inverse (fallback)
#
\(A::AbstractLinearOperator, x::Vector) = adjoint(A)*x
