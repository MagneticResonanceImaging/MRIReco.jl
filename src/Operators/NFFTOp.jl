export NFFTOp
import Base.adjoint

mutable struct NFFTOp{T,F1,F2} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: F1
  ctprod :: F2
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  plan 
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
function NFFTOp(shape::Tuple, tr::Trajectory; nodes=nothing, kargs...)
  nodes==nothing ? nodes=kspaceNodes(tr) : nothing
  plan = NFFTPlan(nodes, shape, 3, 1.25, precompute = NFFT.FULL)

  function produ(x::Vector{T}) where T<:Union{Real,Complex}
    y = nfft(plan,reshape(x[:],shape))
    return vec(y)
  end

  function ctprodu(y::Vector{T}) where T<:Union{Real,Complex}
    x = nfft_adjoint(plan, y[:])
    return vec(x)
  end

  return NFFTOp{ComplexF64,Nothing,Function}(size(nodes,2), prod(shape), false, false
            , produ
            , nothing
            , ctprodu, 0, 0, 0, plan)
end

function adjoint(op::NFFTOp{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod, nothing, op.prod)
end


### Normal Matrix Code ###


struct NFFTNormalOp{S,D} 
  parent::S
  weights::D
  fftplan
  ifftplan
  λ
end

function NFFTNormalOp(S::NFFTOp, W)

  nodes = S.x
  shape = S.N

  # this will not work here since W is an operator and should be just the weights
  λ, ft, ift = diagonalizeOp(S.plan, W)

  return NFFTNormalOp(S, W, ft, ift, λ)
end

function SparsityOperators.normalOperator(S::NFFTOp, W=I)
  return NFFTNormalOp(S,W)
end


function Base.:*(N::NFFTNormalOp, x::AbstractVector)
  shape = N.plan.N

  xL = zeros(T,Tuple(2*collect(shape)))
  xL[1:shape[1],1:shape[2]] = x
  λ = reshape(N.λ,Tuple(2*collect(shape)))
  
  # return vec( (ifftplan*( λ.*(fftplan*xL)))[1:shape[1],1:shape[2]] )
  y = (N.ifftplan*( λ.*(N.fftplan*xL)))[1:shape[1],1:shape[2]]
  
  return vec(y)
end


function diagonalizeOp(p::NFFTPlan, weights=nothing)
  shape = p.N
  nodes = p.x

  # plan the FFTs
  fftplan = plan_fft(zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE)
  ifftplan = plan_ifft(zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE)

  # calculate first column of block Toeplitz matrix
  firstCol = [getMatrixElement(j,1,shape,nodes,weights=weights) for j=1:prod(shape)]
  firstCol = reshape(firstCol, shape)

  # TODO: more efficient implementation using NFFTs
  # e1 = zeros(ComplexF64,shape)
  # e1[1] = 1.
  # firstCol = nfft_adjoint(p, nfft(p,e1)[:].*weights)
  # firstCol = reshape(firstCol, shape)

  # calculate first rows of the leftmost Toeplitz blocks
  firstRow = zeros(ComplexF64,shape[2],shape[1])
  for i=1:shape[2]
     firstRow[i,:] = [ getMatrixElement((i-1)*shape[1]+1,k,shape,nodes,weights=weights) for k=1:shape[1] ] # first row of the relevant blocks
  end

  # firstRow = zeros(ComplexF64,shape[2],shape[1])
  #   for i=1:shape[2]
  #     e1 = zeros(ComplexF64,shape)
  #     e1[1,i] = 1.
  #     firstRow[i,:] = conj( nfft_adjoint(p, nfft(p,e1)[:].*weights)[1:shape[1]] )
  #   end

  # construct FT matrix of the eigentvalues
  eigMat = zeros(ComplexF64, Tuple(2*collect(shape)))
  eigMat[1:shape[1], 1:shape[2]] = firstCol
  eigMat[shape[1]+2:end, 1:shape[2]] = copy(transpose(firstRow[:, end:-1:2]))
  eigMat[1:shape[1], shape[2]+2:end] = firstRow[end:-1:2,:]'
  eigMat[shape[1]+2:end, shape[2]+2:end] =  conj(firstCol[end:-1:2,end:-1:2])

  return vec(fftplan*eigMat), fftplan, ifftplan
end




#
# calculate the matrix element A_{j,k} explicitely
#
function getMatrixElement(j::Int, k::Int, shape::Tuple, nodes::Matrix; weights=nothing)
  elem=0.
  x = mod(k-1,shape[1])-mod(j-1,shape[1])
  y = div(k-1,shape[1])-div(j-1,shape[1])
  if weights != nothing
    for i=1:size(nodes,2)
      elem += exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )*weights[i]
    end
  else
    for i=1:size(nodes,2)
      elem += exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )
    end
  end

  return elem
end