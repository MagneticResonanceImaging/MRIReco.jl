module MRIOperators

import Base: hcat, vcat, \
export hcat, vcat, \, diagOp

using Reexport
using MRIBase
@reexport using SparsityOperators
using LinearAlgebra
using NFFT
using FLoops

import LowRankApprox.psvd
using Distributions
using StatsBase

include("Shutter.jl")
include("Composition.jl")
include("NFFTOp.jl")
include("ExplicitOp.jl")
include("SensitivityOp.jl")
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

function diagOpProd(y::AbstractVector{T}, x::AbstractVector{T}, nrow::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  @floop for i=1:length(ops)
    mul!(view(y,yIdx[i]:yIdx[i+1]-1), ops[i], view(x,xIdx[i]:xIdx[i+1]-1))
  end
  return y
end

function diagOpTProd(y::AbstractVector{T}, x::AbstractVector{T}, ncol::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  @floop for i=1:length(ops)
    mul!(view(y,yIdx[i]:yIdx[i+1]-1), transpose(ops[i]), view(x,xIdx[i]:xIdx[i+1]-1))
  end
  return y
end

function diagOpCTProd(y::AbstractVector{T}, x::AbstractVector{T}, ncol::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  @floop for i=1:length(ops)
    mul!(view(y,yIdx[i]:yIdx[i+1]-1), adjoint(ops[i]), view(x,xIdx[i]:xIdx[i+1]-1))
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
                     (res,x) -> (diagOpProd(res,x,nrow,xIdx,yIdx,ops...)),
                     (res,y) -> (diagOpTProd(res,y,ncol,yIdx,xIdx,ops...)),
                     (res,y) -> (diagOpCTProd(res,y,ncol,yIdx,xIdx,ops...)), 
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
                    (res,x) -> (diagOpProd(res,x,nrow,xIdx,yIdx,ops...)),
                    (res,y) -> (diagOpTProd(res,y,ncol,yIdx,xIdx,ops...)),
                    (res,y) -> (diagOpCTProd(res,y,ncol,yIdx,xIdx,ops...)), 
                     0, 0, 0, false, false, false, S[], S[],
                     ops, true, xIdx, yIdx )

  return Op
end

### Normal Matrix Code ###

mutable struct DiagNormalOp{T,V} <: AbstractLinearOperator{T}
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
  normalOps::V
  idx::Vector{Int64}
  y::Vector{T}
end

LinearOperators.storage_type(op::DiagNormalOp) = typeof(op.Mv5)

function DiagNormalOp(normalOps, N, idx, y::Vector{T}) where {T}
  
  function produ!(y, normalOps, idx, x)
    @floop for i=1:length(normalOps)
       mul!(view(y,idx[i]:idx[i+1]-1), normalOps[i], view(x,idx[i]:idx[i+1]-1))
    end
    return y
  end  

  return DiagNormalOp(N, N, false, false
         , (res,x) -> produ!(res, normalOps, idx, x)
         , nothing
         , nothing
         , 0, 0, 0, false, false, false, T[], T[]
         , normalOps, idx, y)
end

function SparsityOperators.normalOperator(S::DiagOp, W=opEye(eltype(S), size(S,1)))
  weights = W*ones(S.nrow)

  T = promote_type(eltype(S), eltype(W))

  if S.equalOps
    # this optimization is only allowed if all ops are the same

    # we promote the weights to be of the same type as T, which will be required
    # when creating the temporary vector in normalOperator in a later stage
    opInner = normalOperator(S.ops[1], WeightingOp(T.(weights[S.yIdx[1]:S.yIdx[2]-1].^2)))
    op = DiagNormalOp([copy(opInner) for i=1:length(S.ops)], size(S,2), S.xIdx, zeros(T, S.ncol) )
  else
    op = DiagNormalOp([normalOperator(S.ops[i], WeightingOp(T.(weights[S.yIdx[i]:S.yIdx[i+1]-1].^2)))
                     for i in 1:length(S.ops)], size(S,2), S.xIdx, zeros(T, S.ncol) )
  end

  return op
end


end # module