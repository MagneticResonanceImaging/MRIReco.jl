module MRIOperators

using Base: hcat, vcat, \

using Adapt
using Reexport
using MRIBase
@reexport using LinearOperatorCollection
using LinearAlgebra
using NFFT
using NFFT.FFTW
using Wavelets
using FLoops

using StatsBase

include("Shutter.jl")
include("ExplicitOp.jl")
include("SensitivityOp.jl")
include("MapSliceOp.jl")
include("FieldmapNFFTOp.jl")
include("EncodingOp.jl")
include("SparseOp.jl")
include("SubspaceOp.jl")

"""
    hcat(A::AbstractLinearOperator, n::Int)

horizontally concatenates a linear operator `A` n times
"""
function LinearOperators.hcat(A::AbstractLinearOperator, n::Int)
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
function LinearOperators.vcat(A::AbstractLinearOperator, n::Int)
  op = A
  for i = 2:n
    op = LinearOperators.vcat(op, A)
  end
  return op
end

function LinearOperatorCollection.diagOpProd(y::Vector{T}, x::Vector{T}, nrow::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  @floop for i=1:length(ops)
    mul!(view(y,yIdx[i]:yIdx[i+1]-1), ops[i], view(x,xIdx[i]:xIdx[i+1]-1))
  end
  return y
end

function LinearOperatorCollection.diagOpTProd(y::Vector{T}, x::Vector{T}, ncol::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  @floop for i=1:length(ops)
    mul!(view(y,yIdx[i]:yIdx[i+1]-1), transpose(ops[i]), view(x,xIdx[i]:xIdx[i+1]-1))
  end
  return y
end

function LinearOperatorCollection.diagOpCTProd(y::Vector{T}, x::Vector{T}, ncol::Int, xIdx, yIdx, ops :: AbstractLinearOperator...) where T
  @floop for i=1:length(ops)
    mul!(view(y,yIdx[i]:yIdx[i+1]-1), adjoint(ops[i]), view(x,xIdx[i]:xIdx[i+1]-1))
  end
  return y
end

function LinearOperatorCollection.diagNormOpProd!(y::Vector{T}, normalOps, idx, x::Vector{T}) where T
  @floop for i=1:length(normalOps)
    mul!(view(y,idx[i]:idx[i+1]-1), normalOps[i], view(x,idx[i]:idx[i+1]-1))
 end
 return y
end

function Base.copy(S::LinearOperator{T}) where T
  deepcopy(S)
end

fftParams(::Type{<:AbstractArray}) = (;:flags => FFTW.MEASURE)

# https://github.com/JuliaLang/julia/issues/35543
stripParameters(arrT::Type{<:AbstractArray}) = Base.typename(arrT).wrapper

end # module