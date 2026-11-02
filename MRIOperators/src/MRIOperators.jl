module MRIOperators

using Base: hcat, vcat, \

import Base: ∘, parent, getproperty, setproperty!

using Adapt
using Reexport
using MRIBase
@reexport using LinearOperatorCollection
using LinearOperatorCollection.LinearOperators
using LinearAlgebra
using NFFT
using NFFT.AbstractNFFTs
using NFFT.FFTW
using Wavelets
using OhMyThreads

using TSVD: tsvd
using StatsBase

export AbstractMRIOperator
abstract type AbstractMRIOperator{T} <: AbstractLinearOperatorFromCollection{T} end

abstract type WrappedMRIOperator{T, O} <: AbstractMRIOperator{T} end
parent(op::WrappedMRIOperator) = throw(error("$(typeof(op)) must implement `parent`"))
getproperty(op::WrappedMRIOperator, f::Symbol) = getproperty(parent(op), f)
setproperty!(op::WrappedMRIOperator, f::Symbol, v) = setproperty!(parent(op), f, v)
LinearOperators.storage_type(op::WrappedMRIOperator) = LinearOperators.storage_type(parent(op))

# Adapts an input array that a wrapper retains for reconstruction, returning a plain array of the target's
# array type. Adapt may return wrapped arrays (views, structured matrices, ...) or leave an input
# unconverted (array adaptors keep isbits arrays such as StaticArrays' SMatrix). Such results are copied into
# a plain array on the backend that holds them, so device data is copied on the device, and then adapted.
function adaptRetained(to, x::AbstractArray, name)
  target = typeof(adapt(to, Array{eltype(x)}(undef, ntuple(_ -> 0, ndims(x))))).name.wrapper
  y = adapt(to, x)
  y isa target && y isa DenseArray && return y
  dense = similar(y, eltype(y), size(y))
  dense .= y
  dense = adapt(to, dense)
  dense isa target && dense isa DenseArray ||
    throw(ArgumentError("cannot adapt: $name of type $(typeof(x)) could not be copied into a $target array on the target"))
  return dense
end

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