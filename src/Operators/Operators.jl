import Base: hcat, vcat, \
export hcat, vcat, \, diagOp, A_mul_B

include("NFFTOp.jl")
include("ExplicitOp.jl")
include("SensitivityOp.jl")
include("SamplingOp.jl")
include("MapSliceOp.jl")
include("FieldmapNFFTOp.jl")
include("EncodingOp.jl")
include("SparseOp.jl")

"""
    hcat(A::AbstractLinearOperator, n::Int)

horizontally concatenates alinear operator `A` n times
"""
function hcat(A::AbstractLinearOperator, n::Int)
  op = A
  for i = 2:n
    op = LinearOperators.hcat(op, A)
  end
  return op
end

"""
    hcat(A::AbstractLinearOperator, n::Int)

vertically concatenates alinear operator `A` n times
"""
function vcat(A::AbstractLinearOperator, n::Int)
  op = A
  for i = 2:n
    op = LinearOperators.vcat(op, A)
  end
  return op
end

function diagOpProd(x::Vector{T}, nrow::Int, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, nrow)
  xIdx=cumsum(vcat(1,[ops[i].ncol for i=1:length(ops)]))
  yIdx=cumsum(vcat(1,[ops[i].nrow for i=1:length(ops)]))
  @sync for i=1:length(ops)
    Threads.@spawn begin
      y[yIdx[i]:yIdx[i+1]-1] = ops[i]*x[xIdx[i]:xIdx[i+1]-1]
      y[yIdx[i]:yIdx[i+1]-1] = ops[i]*x[xIdx[i]:xIdx[i+1]-1]
    end
  end
  return y
end

function diagOpTProd(x::Vector{T}, ncol::Int, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, ncol)
  xIdx=cumsum(vcat(1,[ops[i].nrow for i=1:length(ops)]))
  yIdx=cumsum(vcat(1,[ops[i].ncol for i=1:length(ops)]))
  @sync for i=1:length(ops)
    Threads.@spawn begin  
      y[yIdx[i]:yIdx[i+1]-1] = transpose(ops[i])*x[xIdx[i]:xIdx[i+1]-1]
    end
  end
  return y
end

function diagOpCTProd(x::Vector{T}, ncol::Int, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, ncol)
  xIdx=cumsum(vcat(1,[ops[i].nrow for i=1:length(ops)]))
  yIdx=cumsum(vcat(1,[ops[i].ncol for i=1:length(ops)]))
  @sync for i=1:length(ops)
    Threads.@spawn begin    
      y[yIdx[i]:yIdx[i+1]-1] = adjoint(ops[i])*x[xIdx[i]:xIdx[i+1]-1]
    end
  end
  return y
end

"""
    diagOp(ops :: AbstractLinearOperator...)

create a bloc-diagonal operator out of the `LinearOperator`s contained in ops
"""
function diagOp(ops :: AbstractLinearOperator...)
  nrow=0
  ncol=0
  S = eltype(ops[1])
  for i = 1:length(ops)
    nrow += ops[i].nrow
    ncol += ops[i].ncol
    S = promote_type(S, eltype(ops[i]))
  end

  return LinearOperator{S}( nrow, ncol, false, false
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
