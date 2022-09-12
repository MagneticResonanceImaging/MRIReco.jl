"""
    `mutable struct CompositeOp{T}`

  struct describing the result of a composition/product of operators.
  Describing the composition using a dedicated type has the advantage
  that the latter can be made copyable. This is particularly relevant for
  multi-threaded code
"""
mutable struct CompositeOp{T} <: AbstractLinearOperator{T}
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
  A
  B
end

"""
    CompositeOp(ops :: AbstractLinearOperator...)

composition/product of two Operators. Differs with * since it can handle normal operator
"""
function CompositeOp(A,B;isWeighting=false)
  nrow=A.nrow
  ncol=B.ncol
  S = eltype(A)

  function produ!(res, x::Vector{T}) where T<:Union{Real,Complex}
    res .= A*(B*x)
  end

  function tprodu!(res, y::Vector{T}) where T<:Union{Real,Complex}
    res .= transposed(B)*(transposed(A)*y)
  end

  function ctprodu!(res, y::Vector{T}) where T<:Union{Real,Complex}
    res .= adjoint(B)*(adjoint(A)*y)
  end

  Op = CompositeOp{S}( nrow, ncol, false, false,
                     produ!,
                     tprodu!,
                     ctprodu!, 
                     0, 0, 0, false, false, false, S[], S[],
                     isWeighting, A, B )

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

struct CompositeNormalOp{S,U} 
  opOuter::S
  normalOpInner::U
end

function SparsityOperators.normalOperator(S::CompositeOp, W=I)
  if S.isWeighting && W==opEye()
    # In this case we are converting the left argument into a 
    # weighting matrix, that is passed to normalOperator
    normalOperator(S.B, S.A)
  else
    return CompositeNormalOp(S.B, normalOperator(S.A, W) )
  end
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
