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
  S = promote_type(eltype(A), eltype(B))
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
function Base.:∘(A::T1, B::T2; isWeighting::Bool=false) where {T1<:AbstractLinearOperator, T2<:AbstractLinearOperator}
  return CompositeOp(A,B;isWeighting=isWeighting)
end

function Base.copy(S::CompositeOp{T}) where T
  A = copy(S.A)
  B = copy(S.B)
  return CompositeOp(A,B; isWeighting=S.isWeighting)
end


### Normal Matrix Code ###
# Left matrix can be build into a normal operator

mutable struct CompositeNormalOp{T,S,U,V} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod! :: Function
  tprod! :: Nothing
  ctprod! :: Nothing
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  args5 :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5 :: Vector{T}
  Mtu5 :: Vector{T}
  opOuter::S
  normalOpInner::U
  tmp::V
end

LinearOperators.storage_type(op::CompositeNormalOp) = typeof(op.Mv5)


function CompositeNormalOp(opOuter, normalOpInner, tmp::Vector{T}) where T

  function produ!(y, opOuter, normalOpInner, tmp, x)
    mul!(tmp, opOuter, x)
    mul!(tmp, normalOpInner, tmp) # This can be dangerous. We might need to create two tmp vectors
    return mul!(y, adjoint(opOuter), tmp)
  end

  return CompositeNormalOp(size(opOuter,2), size(opOuter,2), false, false
         , (res,x) -> produ!(res, opOuter, normalOpInner, tmp, x)
         , nothing
         , nothing
         , 0, 0, 0, false, false, false, T[], T[]
         , opOuter, normalOpInner, tmp)
end


function LinearOperatorCollection.normalOperator(S::CompositeOp, W=opEye(eltype(S),size(S,1)))
  T = promote_type(eltype(S.A), eltype(S.B), eltype(W))
  if S.isWeighting #&& typeof(W) <: opEye
    # In this case we are converting the left argument into a 
    # weighting matrix, that is passed to normalOperator
    normalOperator(S.B, S.A)
  else
    tmp = Vector{T}(undef, size(S.A, 2))
    return CompositeNormalOp(S.B, normalOperator(S.A, W), tmp)
  end
end

function Base.copy(S::CompositeNormalOp) 
  opOuter = copy(S.opOuter)
  opInner = copy(S.normalOpInner)
  tmp = copy(S.tmp)
  return CompositeNormalOp(opOuter, opInner, tmp)
end
