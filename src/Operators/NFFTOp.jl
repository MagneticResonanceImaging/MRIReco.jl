export NFFTOp
import Base.adjoint

mutable struct NFFTOp{T} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: Nothing
  ctprod :: Function
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  plan 
  toeplitz :: Bool
end

"""
    NFFTOp(shape::Tuple, tr::Trajectory; nodes=nothing, kargs...)

generates a `NFFTOp` which evaluates the MRI Fourier signal encoding operator using the NFFT.

# Arguments:
* `shape::NTuple{D,Int64}`  - size of image to encode/reconstruct
* `tr::Trajectory`          - Trajectory with the kspace nodes to sample
* (`nodes=nothing`)         - Array containg the trajectory nodes (redundant)
* (`kargs`)                   - additional keyword arguments
"""
function NFFTOp(shape::Tuple, tr::Trajectory; nodes=nothing, toeplitz=false, 
                oversamplingFactor=1.25, kernelSize=3, kargs...)
  nodes==nothing ? nodes=kspaceNodes(tr) : nothing
  # plan = NFFTPlan(nodes, shape, m=kernelSize, σ=oversamplingFactor, precompute=NFFT.FULL)
  plan = plan_nfft(nodes, shape, m=kernelSize, σ=oversamplingFactor, precompute=NFFT.FULL)

  return NFFTOp{ComplexF64}(size(nodes,2), prod(shape), false, false
            , x->produ(plan,x)
            , nothing
            , y->ctprodu(plan,y), 0, 0, 0, plan, toeplitz)
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
  return NFFTOp{ComplexF64}(size(plan.x,2), prod(plan.N), false, false
         , x->produ(plan,x)
         , nothing
         , y->ctprodu(plan,y), 0, 0, 0, plan, S.toeplitz)
end

### Normal Matrix Code ###

struct NFFTNormalOp{S,D} 
  shape::S
  weights::D
  fftplan
  ifftplan
  λ::Array{ComplexF64}
  xL::Matrix{ComplexF64}
end

function Base.copy(S::NFFTNormalOp)
  fftplan = plan_fft(zeros(ComplexF64, Tuple(2*collect(S.shape)));flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(ComplexF64, Tuple(2*collect(S.shape)));flags=FFTW.MEASURE)
  return NFFTNormalOp(S.shape, S.weights, fftplan, ifftplan, S.λ, S.xL)
end

function NFFTNormalOp(S::NFFTOp, W)
  shape = S.plan.N
  weights = W*ones(S.nrow)

  λ, ft, ift = diagonalizeOp(S.plan, weights)
  xL = zeros(ComplexF64,2*shape[1], 2*shape[2])

  return NFFTNormalOp(shape, W, ft, ift, λ, xL)
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

  # xL = zeros(T,Tuple(2*collect(shape)))
  # xL[1:shape[1],1:shape[2]] = x
  # λ = reshape(S.λ,Tuple(2*collect(shape)))
  
  # y = (S.ifftplan*( λ.*(S.fftplan*xL)))[1:shape[1],1:shape[2]]

  S.xL .= 0
  S.xL[1:shape[1],1:shape[2]] .= reshape(x,shape)
  λ = reshape(S.λ,Tuple(2*collect(shape)))
  
  y = (S.ifftplan*( λ.*(S.fftplan*S.xL)))[1:shape[1],1:shape[2]]
  
  return vec(y)
end


function diagonalizeOp(p::NFFT.NFFTPlan, weights=nothing)
  shape = p.N
  nodes = p.x

  # plan the FFTs
  fftplan = plan_fft(zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE)

  # calculate first column of block Toeplitz matrix
   e1 = zeros(ComplexF64,shape)
   e1[1] = 1.
   firstCol = nfft_adjoint(p, nfft(p,e1)[:].*weights)
   firstCol = reshape(firstCol, shape)

  # calculate first rows of the leftmost Toeplitz blocks
  p2 = plan_nfft([-1,1].*nodes, shape, m=3, σ=1.25)
  firstRow = nfft_adjoint(p2, nfft(p2,e1)[:].*weights)
  firstRow = reshape(firstRow, shape)

  # construct FT matrix of the eigentvalues
  # here a row and column in the center is filled with zeros.
  # this is done to ensure that the FFT can operate
  # on arrays with even size in each dimension
  eigMat = zeros(ComplexF64, Tuple(2*collect(shape)))
  eigMat[1:shape[1], 1:shape[2]] .= firstCol
  eigMat[shape[1]+2:end, 1:shape[2]] .= firstRow[end:-1:2,:]
  eigMat[1:shape[1], shape[2]+2:end] .= conj.(firstRow[:,end:-1:2])
  eigMat[shape[1]+2:end, shape[2]+2:end] .=  conj.(firstCol[end:-1:2,end:-1:2])

  return vec(fftplan*eigMat), fftplan, ifftplan
end




#
# calculate the matrix element A_{j,k} explicitely
#
#function getMatrixElement(j::Int, k::Int, shape::Tuple, nodes::Matrix; weights=nothing)
#  elem=0.
#  x = mod(k-1,shape[1])-mod(j-1,shape[1])
#  y = div(k-1,shape[1])-div(j-1,shape[1])
#  if weights != nothing
#    for i=1:size(nodes,2)
#      elem += exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )*weights[i]
#    end
#  else
#    for i=1:size(nodes,2)
#      elem += exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )
#    end
#  end
#
#  return elem
#end
