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


mutable struct DiagOp{T} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: Function
  ctprod :: Function
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
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
                     x->diagOpProd(x,nrow,xIdx,yIdx,ops...),
                     y->diagOpTProd(y,ncol,yIdx,xIdx,ops...),
                     y->diagOpCTProd(y,ncol,yIdx,xIdx,ops...), 0, 0, 0, 
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
                     x->diagOpProd(x,nrow,xIdx,yIdx,ops...),
                     y->diagOpTProd(y,ncol,yIdx,xIdx,ops...),
                     y->diagOpCTProd(y,ncol,yIdx,xIdx,ops...), 0, 0, 0, 
                     ops, true, xIdx, yIdx )

  return Op
end



### Normal Matrix Code ###

struct DiagNormalOp
  ops
  normalOps
  nrow
  ncol
  idx
  y
end

function SparsityOperators.normalOperator(S::DiagOp, W=I)
  weights = W*ones(S.nrow)

  if S.equalOps
    # this opimization is only allowed if all ops are the same
    opInner = normalOperator(S.ops[1], WeightingOp(weights[S.yIdx[1]:S.yIdx[2]-1].^2))
    op = DiagNormalOp(S.ops, [copy(opInner) for i=1:length(S.ops)], S.ncol, S.ncol, S.xIdx, zeros(eltype(S), S.ncol) )
  else
    op = DiagNormalOp(S.ops, [normalOperator(S.ops[i], WeightingOp(weights[S.yIdx[i]:S.yIdx[i+1]-1].^2)) 
                     for i in 1:length(S.ops)], S.ncol, S.ncol, S.xIdx, zeros(eltype(S), S.ncol) )
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
       y[idx[i]:idx[i+1]-1] = ops[i]*x[idx[i]:idx[i+1]-1]
    end
  end
  return
end

mutable struct ProdOp{T} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: Function
  ctprod :: Function
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  isWeighting :: Bool
  A
  B
end


"""
    prodOp(ops :: AbstractLinearOperator...)

product of two Operators. Differs with * since it can handle normal operator
"""
function prodOp(A,B;isWeighting=false)
  nrow=A.nrow
  ncol=B.ncol
  S = eltype(A)

  function produ(x::Vector{T}) where T<:Union{Real,Complex}
    return A*(B*x)
  end

  function tprodu(y::Vector{T}) where T<:Union{Real,Complex}
    return transposed(B)*(transposed(A)*y)
  end

  function ctprodu(y::Vector{T}) where T<:Union{Real,Complex}
    return adjoint(B)*(adjoint(A)*y)
  end

  Op = ProdOp{S}( nrow, ncol, false, false,
                     produ,
                     tprodu,
                     ctprodu, 0, 0, 0, isWeighting, A, B )

  return Op
end




### Normal Matrix Code ###
# Left matrix can be build into a normal operator

struct ProdNormalOp{S,U} 
  opOuter::S
  normalOpInner::U
end

function SparsityOperators.normalOperator(S::ProdOp, W=I)
  if S.isWeighting && W==I
    normalOperator(S.B, S.A)
  else
    return ProdNormalOp(S.B, normalOperator(S.A, W) )
  end
end

function Base.:*(S::ProdNormalOp, x::AbstractVector)
  return adjoint(S.opOuter)*(S.normalOpInner*(S.opOuter*x))
end


# implement A_mul_B for the product
A_mul_B(A::AbstractLinearOperator, x::Vector) = A*x

#
# use hermitian conjugate as an estimate for the inverse (fallback)
#
\(A::AbstractLinearOperator, x::Vector) = adjoint(A)*x
