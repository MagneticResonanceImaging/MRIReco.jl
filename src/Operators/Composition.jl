"""
    `mutable struct CompositeOp{T}`

  struct describing the result of a composition/product of operators.
  Describing the composition using a dedicated type has the advantage
  that the latter can be made copyable. This is particularly relevant for
  multi-threaded code
"""
mutable struct CompositeOp{T,U,V} <: AbstractLinearOperator{T}
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
  isWeighting :: Bool
  A::U
  B::V
  tmp::Vector{T}
end

"""
    CompositeOp(ops :: AbstractLinearOperator...)

composition/product of two Operators. Differs with * since it can handle normal operator
"""
function CompositeOp(A,B;isWeighting=false)
  nrow = A.nrow
  ncol = B.ncol
  S = eltype(A)
  tmp_ = Vector{S}(undef, B.nrow)

  function produ!(res, x::AbstractVector{T}, tmp) where T<:Union{Real,Complex}
    mul!(tmp, B, x)
    return mul!(res, A, tmp)
  end

  function tprodu!(res, y::AbstractVector{T}, tmp) where T<:Union{Real,Complex}
    mul!(tmp, transpose(A), y)
    return mul!(res, transpose(B), tmp)
  end

  function ctprodu!(res, y::AbstractVector{T}, tmp) where T<:Union{Real,Complex}
    mul!(tmp, adjoint(A), y)
    return mul!(res, adjoint(B), tmp)
  end

  Op = CompositeOp( nrow, ncol, false, false,
                     (res,x) -> produ!(res,x,tmp_),
                     (res,y) -> tprodu!(res,y,tmp_),
                     (res,y) -> ctprodu!(res,y,tmp_), 
                     0, 0, 0, false, false, false, S[], S[],
                     isWeighting, A, B, tmp_)

  return Op
end

LinearOperators.storage_type(op::CompositeOp) = typeof(op.Mv5)

"""
∘(A::T1, B::T2; isWeighting::Bool=false) where {T1<:AbstractLinearOperator, T2<:AbstractLinearOperator}

  composition of two operators.
"""
function ∘(A::T1, B::T2; isWeighting::Bool=false) where {T1<:AbstractLinearOperator, T2<:AbstractLinearOperator}
  return CompositeOp(A,B;isWeighting=isWeighting)
end

function Base.copy(S::CompositeOp{T}) where T
  A = copy(S.A)
  B = copy(S.B)
  return CompositeOp(A,B; isWeighting=S.isWeighting)
end


### Normal Matrix Code ###
# Left matrix can be build into a normal operator

struct CompositeNormalOp{S,U,V} 
  opOuter::S
  normalOpInner::U
  tmp::V
end

function SparsityOperators.normalOperator(S::CompositeOp, W=opEye(eltype(S),size(S,1)))
  if S.isWeighting #&& typeof(W) <: opEye
    # In this case we are converting the left argument into a 
    # weighting matrix, that is passed to normalOperator
    normalOperator(S.B, S.A)
  else
    tmp = Vector{eltype(S.B)}(undef, size(S.A, 2))
    return CompositeNormalOp(S.B, normalOperator(S.A, W), tmp)
  end
end

function LinearAlgebra.mul!(y, S::CompositeNormalOp, x)
  mul!(S.tmp, S.opOuter, x)
  mul!(S.tmp, S.normalOpInner, S.tmp) # This can be dangerous. We might need to create two tmp vectors
  return mul!(y, adjoint(S.opOuter), S.tmp)
end

function Base.:*(S::CompositeNormalOp, x::AbstractVector)
  return adjoint(S.opOuter)*(S.normalOpInner*(S.opOuter*x))
end


# implement A_mul_B for the product
A_mul_B(A::AbstractLinearOperator, x::Vector) = A*x

function Base.copy(S::CompositeNormalOp{T,U}) where {T,U}
  opOuter = copy(S.opOuter)
  opInner = copy(S.normalOpInner)
  return CompositeNormalOp(opOuter, opInner)
end
