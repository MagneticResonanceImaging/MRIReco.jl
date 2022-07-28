export NFFTOp
import Base.adjoint

mutable struct NFFTOp{T} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod! :: Function
  tprod! :: Nothing
  ctprod! :: Function
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  args5 :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5 :: Vector{T}
  Mtu5 :: Vector{T}
  plan
  toeplitz :: Bool
end

# LinearOperators.storage_type(op::NFFTOp) = typeof(op.Mv5)

"""
    NFFTOp(shape::Tuple, tr::Trajectory; kargs...)
    NFFTOp(shape::Tuple, tr::AbstractMatrix; kargs...)

generates a `NFFTOp` which evaluates the MRI Fourier signal encoding operator using the NFFT.

# Arguments:
* `shape::NTuple{D,Int64}`  - size of image to encode/reconstruct
* `tr`                      - Either a `Trajectory` object, or a `ND x Nsamples` matrix for an ND-dimenensional (e.g. 2D or 3D) NFFT with `Nsamples` k-space samples
* (`nodes=nothing`)         - Array containg the trajectory nodes (redundant)
* (`kargs`)                 - additional keyword arguments
"""
function NFFTOp(shape::Tuple, tr::AbstractMatrix{T}; toeplitz=false, oversamplingFactor=1.25, kernelSize=3, kargs...) where {T}

  plan = plan_nfft(tr, shape, m=kernelSize, σ=oversamplingFactor, precompute=NFFT.FULL)

  return NFFTOp{Complex{T}}(size(tr,2), prod(shape), false, false
            , (res,x) -> (res .= produ(plan,x))
            , nothing
            , (res,y) -> (res .= ctprodu(plan,y))
            , 0, 0, 0, false, false, false, Complex{T}[], Complex{T}[]
            , plan, toeplitz)
end


function NFFTOp(shape::Tuple, tr::Trajectory; toeplitz=false, oversamplingFactor=1.25, kernelSize=3, kargs...)
  return NFFTOp(shape, kspaceNodes(tr); toeplitz=toeplitz, oversamplingFactor=oversamplingFactor, kernelSize=kernelSize, kargs...)
end

function produ(plan::NFFT.NFFTPlan, x::Vector{T}) where T<:Union{Real,Complex}
  y = plan * reshape(x[:],plan.N)
  return vec(y)
end

function ctprodu(plan::NFFT.NFFTPlan, y::Vector{T}) where T<:Union{Real,Complex}
  x = adjoint(plan) * y[:]
  return vec(x)
end


function Base.copy(S::NFFTOp{T}) where {T}
  plan = copy(S.plan)
  return NFFTOp{T}(size(k,2), prod(plan.N), false, false
              , (res,x) -> (res .= produ(plan,x))
              , nothing
              , (res,y) -> (res .= ctprodu(plan,y))
              , 0, 0, 0, false, false, false, T[], T[]
              , plan, S.toeplitz)
end



#########################################################################
### Toeplitz Operator ###
#########################################################################
struct NFFTToeplitzNormalOp{T,D,W}
  shape::NTuple{D,Int}
  weights::W
  fftplan
  ifftplan
  λ::Array{T}
  xL1::Array{T,D}
  xL2::Array{T,D}
end


function NFFTToeplitzNormalOp(S::NFFTOp{T}, W::UniformScaling=I) where {T}
  shape = S.plan.N

  # plan the FFTs
  fftplan  = plan_fft( zeros(T, 2 .* shape);flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(T, 2 .* shape);flags=FFTW.MEASURE)

  λ = calculateToeplitzKernel(shape, S.k; m = S.plan.params.m, σ = S.plan.params.σ, window = S.plan.params.window, LUTSize = S.plan.params.LUTSize, fftplan = fftplan)

  xL1 = Array{T}(undef, 2 .* shape)
  xL2 = similar(xL1)

  return NFFTToeplitzNormalOp(shape, W, fftplan, ifftplan, λ, xL1, xL2)
end

function SparsityOperators.normalOperator(S::NFFTOp, W::UniformScaling=I)
  if S.toeplitz
    return NFFTToeplitzNormalOp(S,W)
  else
    return NormalOp(S,W)
  end
end

function SparsityOperators.normalOperator(S::NFFTOp, W)
  if S.toeplitz
    @warn "Topelitz with non-Uniform scaling is currently not implemented. Using a non-Toeplitz implemenation instead."
  end
  return NormalOp(S,W)
end


function LinearAlgebra.mul!(y, S::NFFTToeplitzNormalOp, b)
  S.xL1 .= 0
  b = reshape(b, S.shape)

  S.xL1[CartesianIndices(b)] .= b
  mul!(S.xL2, S.fftplan, S.xL1)
  S.xL2 .*= S.λ
  mul!(S.xL1, S.ifftplan, S.xL2)

  y .= vec(S.xL1[CartesianIndices(b)])
  return y
end


Base.:*(S::NFFTToeplitzNormalOp, b::AbstractVector) = mul!(similar(b), S, b)
Base.size(S::NFFTToeplitzNormalOp) = S.shape
Base.size(S::NFFTToeplitzNormalOp, dim) = S.shape[dim]
Base.eltype(::Type{NFFTToeplitzNormalOp{T,D,W}}) where {T,D,W} = T

function Base.copy(A::NFFTToeplitzNormalOp{T,D,W}) where {T,D,W}
  fftplan  = plan_fft( zeros(T, 2 .* A.shape); flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(T, 2 .* A.shape); flags=FFTW.MEASURE)
  return NFFTToeplitzNormalOp(A.shape, A.weights, fftplan, ifftplan, A.λ, A.xL1, A.xL2)
end
