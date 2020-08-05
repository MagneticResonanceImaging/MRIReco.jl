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
  plan = NFFTPlan(nodes, shape, 3, 1.25, precompute=NFFT.FULL)

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
            , ctprodu, 0, 0, 0)
end

function adjoint(op::NFFTOp{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod, nothing, op.prod)
end
