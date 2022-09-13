import Base: hcat, vcat, \
export hcat, vcat, \, diagOp, A_mul_B

include("Composition.jl")
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

function diagOpProd(x::Vector{T}, nrow::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, nrow)
  @sync for i=1:length(ops)
    Threads.@spawn begin
      y[yIdx[i]:yIdx[i+1]-1] = ops[i]*x[xIdx[i]:xIdx[i+1]-1]
    end
  end
  return y
end

function diagOpTProd(x::Vector{T}, ncol::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, ncol)
  @sync for i=1:length(ops)
    Threads.@spawn begin  
      y[yIdx[i]:yIdx[i+1]-1] = transpose(ops[i])*x[xIdx[i]:xIdx[i+1]-1]
    end
  end
  return y
end

function diagOpCTProd(x::Vector{T}, ncol::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  y = Vector{T}(undef, ncol)
  @sync for i=1:length(ops)
    Threads.@spawn begin    
      y[yIdx[i]:yIdx[i+1]-1] = adjoint(ops[i])*x[xIdx[i]:xIdx[i+1]-1]
    end
  end
  return y
end

function Base.copy(S::LinearOperator{T}) where T
  deepcopy(S)
end

mutable struct DiagOp{T} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod! :: Function
  tprod! :: Function
  ctprod! :: Function
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  args5 :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5 :: Vector{T}
  Mtu5 :: Vector{T}
  ops
  equalOps :: Bool
  xIdx :: Vector{Int}
  yIdx :: Vector{Int}
end


"""
    diagOp(ops :: AbstractLinearOperator...)

create a bloc-diagonal operator out of the `LinearOperator`s contained in ops
"""
function diagOp(ops :: AbstractLinearOperator...)
  nrow = 0
  ncol = 0
  S = eltype(ops[1])
  for i = 1:length(ops)
    nrow += ops[i].nrow
    ncol += ops[i].ncol
    S = promote_type(S, eltype(ops[i]))
  end

  xIdx = cumsum(vcat(1,[ops[i].ncol for i=1:length(ops)]))
  yIdx = cumsum(vcat(1,[ops[i].nrow for i=1:length(ops)]))

  Op = DiagOp{S}( nrow, ncol, false, false,
                     (res,x) -> (res .= diagOpProd(x,nrow,xIdx,yIdx,ops...)),
                     (res,y) -> (res .= diagOpTProd(y,ncol,yIdx,xIdx,ops...)),
                     (res,y) -> (res .= diagOpCTProd(y,ncol,yIdx,xIdx,ops...)), 
                     0, 0, 0, false, false, false, S[], S[],
                     [ops...], false, xIdx, yIdx)

  return Op
end

function diagOp(op::AbstractLinearOperator, N=1)
  nrow = N*op.nrow
  ncol = N*op.ncol
  S = eltype(op)
  ops = [copy(op) for n=1:N]

  xIdx = cumsum(vcat(1,[ops[i].ncol for i=1:length(ops)]))
  yIdx = cumsum(vcat(1,[ops[i].nrow for i=1:length(ops)]))

  Op = DiagOp{S}( nrow, ncol, false, false,
                    (res,x) -> (res .= diagOpProd(x,nrow,xIdx,yIdx,ops...)),
                    (res,y) -> (res .= diagOpTProd(y,ncol,yIdx,xIdx,ops...)),
                    (res,y) -> (res .= diagOpCTProd(y,ncol,yIdx,xIdx,ops...)), 
                     0, 0, 0, false, false, false, S[], S[],
                     ops, true, xIdx, yIdx )

  return Op
end



### Normal Matrix Code ###

struct DiagNormalOp{U,V,T}
  ops::U
  normalOps::V
  nrow::Int64
  ncol::Int64
  idx::Vector{Int64}
  y::Vector{T}
end

function SparsityOperators.normalOperator(S::DiagOp, W=I)
  weights = W*ones(S.nrow)

  T = eltype(S)

  if S.equalOps
    # this opimization is only allowed if all ops are the same
    opInner = normalOperator(S.ops[1], WeightingOp(weights[S.yIdx[1]:S.yIdx[2]-1].^2))
    op = DiagNormalOp(S.ops, [copy(opInner) for i=1:length(S.ops)], S.ncol, S.ncol, S.xIdx, zeros(T, S.ncol) )
  else
    op = DiagNormalOp(S.ops, [normalOperator(S.ops[i], WeightingOp(weights[S.yIdx[i]:S.yIdx[i+1]-1].^2)) 
                     for i in 1:length(S.ops)], S.ncol, S.ncol, S.xIdx, zeros(T, S.ncol) )
  end

  return op
end

function Base.:*(S::DiagNormalOp, x::AbstractVector{T}) where T
  _produ_diagnormalop(S.normalOps, S.idx, x, S.y) 
  return S.y
end

function _produ_diagnormalop(ops, idx, x, y)
  @sync for i=1:length(ops)
    Threads.@spawn begin
       mul!(view(y,idx[i]:idx[i+1]-1), ops[i], view(x,idx[i]:idx[i+1]-1))
       #y[idx[i]:idx[i+1]-1] = ops[i]*x[idx[i]:idx[i+1]-1]
    end
  end
  return
end

#
# use hermitian conjugate as an estimate for the inverse (fallback)
#
\(A::AbstractLinearOperator, x::Vector) = adjoint(A)*x
