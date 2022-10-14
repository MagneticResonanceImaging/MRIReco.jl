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

LinearOperators.storage_type(op::NFFTOp) = typeof(op.Mv5)

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

  plan = plan_nfft(tr, shape, m=kernelSize, σ=oversamplingFactor, precompute=NFFT.TENSOR,
		                          fftflags=FFTW.ESTIMATE, blocking=true)

  return NFFTOp{Complex{T}}(size(tr,2), prod(shape), false, false
            , (res,x) -> produ!(res,plan,x)
            , nothing
            , (res,y) -> ctprodu!(res,plan,y)
            , 0, 0, 0, false, false, false, Complex{T}[], Complex{T}[]
            , plan, toeplitz)
end


function NFFTOp(shape::Tuple, tr::Trajectory; toeplitz=false, oversamplingFactor=1.25, kernelSize=3, kargs...)
  return NFFTOp(shape, kspaceNodes(tr); toeplitz=toeplitz, oversamplingFactor=oversamplingFactor, kernelSize=kernelSize, kargs...)
end

function produ!(y::AbstractVector, plan::NFFT.NFFTPlan, x::AbstractVector) 
  mul!(y, plan, reshape(x,plan.N))
end

function ctprodu!(x::AbstractVector, plan::NFFT.NFFTPlan, y::AbstractVector)
  mul!(reshape(x, plan.N), adjoint(plan), y)
end


function Base.copy(S::NFFTOp{T}) where {T}
  plan = copy(S.plan)
  return NFFTOp{T}(size(plan.k,2), prod(plan.N), false, false
              , (res,x) -> produ!(res,plan,x)
              , nothing
              , (res,y) -> ctprodu!(res,plan,y)
              , 0, 0, 0, false, false, false, T[], T[]
              , plan, S.toeplitz)
end



#########################################################################
### Toeplitz Operator ###
#########################################################################

mutable struct NFFTToeplitzNormalOp{T,D,W} <: AbstractLinearOperator{T}
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
  shape::NTuple{D,Int}
  weights::W
  fftplan
  ifftplan
  λ::Array{T}
  xL1::Array{T,D}
  xL2::Array{T,D}
end

LinearOperators.storage_type(op::NFFTToeplitzNormalOp) = typeof(op.Mv5)

function NFFTToeplitzNormalOp(shape, W, fftplan, ifftplan, λ, xL1::Array{T,D}, xL2::Array{T,D}) where {T,D}

  function produ!(y, shape, fftplan, ifftplan, λ, xL1, xL2, x)
    xL1 .= 0
    x = reshape(x, shape)
  
    xL1[CartesianIndices(x)] .= x
    mul!(xL2, fftplan, xL1)
    xL2 .*= λ
    mul!(xL1, ifftplan, xL2)
  
    y .= vec(xL1[CartesianIndices(x)])
    return y
  end

  return NFFTToeplitzNormalOp(prod(shape), prod(shape), false, false
         , (res,x) -> produ!(res, shape, fftplan, ifftplan, λ, xL1, xL2, x)
         , nothing
         , nothing
         , 0, 0, 0, false, false, false, T[], T[]
         , shape, W, fftplan, ifftplan, λ, xL1, xL2)
end

function NFFTToeplitzNormalOp(S::NFFTOp{T}, W=opEye(T,size(S,1))) where {T}
  shape = S.plan.N

  # plan the FFTs
  fftplan  = plan_fft( zeros(T, 2 .* shape);flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(T, 2 .* shape);flags=FFTW.MEASURE)

  # TODO extend the following function by weights
  # λ = calculateToeplitzKernel(shape, S.plan.k; m = S.plan.params.m, σ = S.plan.params.σ, window = S.plan.params.window, LUTSize = S.plan.params.LUTSize, fftplan = fftplan)

  shape_os = 2 .* shape
  p = plan_nfft(typeof(S.plan.k), S.plan.k, shape_os; m = S.plan.params.m, σ = S.plan.params.σ,
		precompute=NFFT.POLYNOMIAL, fftflags=FFTW.ESTIMATE, blocking=true)
  eigMat = adjoint(p) * ( W  * ones(T, size(S.plan.k,2)))
  λ = fftplan * fftshift(eigMat)

  xL1 = Array{T}(undef, 2 .* shape)
  xL2 = similar(xL1)

  return NFFTToeplitzNormalOp(shape, W, fftplan, ifftplan, λ, xL1, xL2)
end

function SparsityOperators.normalOperator(S::NFFTOp{T}, W=opEye(T,size(S,1))) where T
  if S.toeplitz
    return NFFTToeplitzNormalOp(S,W)
  else
    return NormalOp(S,W)
  end
end

function Base.copy(A::NFFTToeplitzNormalOp{T,D,W}) where {T,D,W}
  fftplan  = plan_fft( zeros(T, 2 .* A.shape); flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(T, 2 .* A.shape); flags=FFTW.MEASURE)
  return NFFTToeplitzNormalOp(A.shape, A.weights, fftplan, ifftplan, A.λ, copy(A.xL1), copy(A.xL2))
end
