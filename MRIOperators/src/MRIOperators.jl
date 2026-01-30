module MRIOperators

using Base: hcat, vcat, \

using Adapt
using Reexport
using MRIBase
@reexport using LinearOperatorCollection
using LinearAlgebra
using NFFT
using NFFT.AbstractNFFTs
using NFFT.FFTW
using Wavelets
using OhMyThreads

using TSVD: tsvd
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

function Base.copy(S::LinearOperator{T}) where T
  deepcopy(S)
end

fftParams(T::Type{<:AbstractArray}) = (;:flags => FFTW.MEASURE)
nfftParams(T::Type{<:AbstractArray}) = nfftParams(AbstractNFFTs.active_backend(), T)
nfftParams(::AbstractNFFTBackend, ::Type{<:AbstractArray}) = (;:fftflags => FFTW.MEASURE)


# https://github.com/JuliaLang/julia/issues/35543
stripParameters(arrT::Type{<:AbstractArray}) = Base.typename(arrT).wrapper

end # module