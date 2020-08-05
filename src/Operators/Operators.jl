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

  Op = DiagOp{S}( nrow, ncol, false, false,
                     x->diagOpProd(x,nrow,ops...),
                     y->diagOpTProd(y,ncol,ops...),
                     y->diagOpCTProd(y,ncol,ops...), 0, 0, 0, [ops...] )

  return Op
end



### Normal Matrix Code ###

struct DiagNormalOp
  ops
  normalOps
  nrow
  ncol
end

function SparsityOperators.normalOperator(S::DiagOp, W=I)
  weights = W*ones(S.nrow)
  yIdx = cumsum(vcat(1,[S.ops[i].nrow for i=1:length(S.ops)]))

  # this line is extremly expensive -> redundant work
  #@time op = DiagNormalOp(S.ops, [normalOperator(S.ops[i], WeightingOp(weights[yIdx[i]:yIdx[i+1]-1].^2)) 
  #                   for i in 1:length(S.ops)], S.ncol, S.ncol )

  # this opimization is only allow if all ops are the same
  t = @elapsed opInner = normalOperator(S.ops[1], WeightingOp(weights[yIdx[1]:yIdx[2]-1].^2))
  @info "Time so build normalOp: $t seconds"
  op = DiagNormalOp(S.ops, [copy(opInner) for i=1:length(S.ops)], S.ncol, S.ncol )

  return op
end

function Base.:*(S::DiagNormalOp, x::AbstractVector{T}) where T
  y = Vector{T}(undef, S.nrow)
  xIdx=cumsum(vcat(1,[S.ops[i].ncol for i=1:length(S.ops)]))
  yIdx=cumsum(vcat(1,[S.ops[i].ncol for i=1:length(S.ops)]))
  @sync for i=1:length(S.ops)
    Threads.@spawn begin
      y[yIdx[i]:yIdx[i+1]-1] = S.normalOps[i]*x[xIdx[i]:xIdx[i+1]-1]
    end
  end
  return y  
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
