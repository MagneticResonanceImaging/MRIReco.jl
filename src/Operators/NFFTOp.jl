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
  y = nfft(plan,reshape(x[:],plan.N))
  return vec(y)
end

function ctprodu(plan::NFFT.NFFTPlan, y::Vector{T}) where T<:Union{Real,Complex}
  x = nfft_adjoint(plan, y[:])
  return vec(x)
end


function Base.copy(S::NFFTOp)
  plan = copy(S.plan)
  return NFFTOp{Complex{T}}(size(plan.x,2), prod(plan.N), false, false
              , (res,x) -> (res .= produ(plan,x))
              , nothing
              , (res,y) -> (res .= ctprodu(plan,y))
              , 0, 0, 0, false, false, false, Complex{T}[], Complex{T}[]
              , plan, S.toeplitz)
end

### Normal Matrix Code ###

struct NFFTNormalOp{T,D,W}
  shape::NTuple{D,Int}
  weights::W
  fftplan
  ifftplan
  λ::Array{T}
  xL::Array{T,D}
end

function Base.copy(A::NFFTNormalOp{T,D,W}) where {T,D,W}
  fftplan  = plan_fft( zeros(T, 2 .* A.shape); flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(T, 2 .* A.shape); flags=FFTW.MEASURE)
  return NFFTNormalOp(A.shape, A.weights, fftplan, ifftplan, A.λ, A.xL)
end

function Base.size(S::NFFTNormalOp)
  return S.shape
end

function Base.size(S::NFFTNormalOp, dim)
  return S.shape[dim]
end

function LinearAlgebra.mul!(x, S::NFFTNormalOp, b)
  x .= S * b
  return x
end


function NFFTNormalOp(S::NFFTOp{T}, W) where {T}
  shape = S.plan.N
  weights = W * ones(T, S.nrow)


  # plan the FFTs
  fftplan  = plan_fft( zeros(T, 2 .* shape);flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(T, 2 .* shape);flags=FFTW.MEASURE)

  λ = calculateToeplitzKernel(shape, S.plan.x; m = S.plan.params.m, σ = S.plan.params.σ, window = S.plan.params.window, LUTSize = S.plan.params.LUTSize, fftplan = fftplan)

  xL = zeros(T, 2 .* shape)

  return NFFTNormalOp(shape, W, fftplan, ifftplan, λ, xL)
end

function SparsityOperators.normalOperator(S::NFFTOp, W=I)
  if S.toeplitz
    return NFFTNormalOp(S,W)
  else
    return NormalOp(S,W)
  end
end


function Base.:*(S::NFFTNormalOp, x::AbstractVector{T}) where T
  shape = S.shape

  S.xL .= 0
  x = reshape(x,shape)

  S.xL[CartesianIndices(x)] .= x

  λ = reshape(S.λ,Tuple(2*collect(shape)))

  y = (S.ifftplan*( λ.*(S.fftplan*S.xL)))[CartesianIndices(x)]

  return vec(y)
end