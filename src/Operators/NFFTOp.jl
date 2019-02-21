export NFFTOp
import Base.adjoint
import LinearOperators.FuncOrNothing

mutable struct NFFTOp{T,F1<:FuncOrNothing,F2<:FuncOrNothing,F3<:FuncOrNothing} <:
                      AbstractLinearOperator{T,F1,F2,F3}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: F1
  ctprod :: F2
  inv :: F3
  density::Vector{Float64}
  circShutter::Bool
end

#
# Linear Operator to perform NFFT
#
function NFFTOp(shape::Tuple, tr::AbstractTrajectory; nodes=nothing, symmetrize=true)
  nodes==nothing ? nodes=kspaceNodes(tr) : nothing
  plan = NFFTPlan(nodes, shape, 3, 1.25)
  density = sdc(plan)

  function produ(x::Vector{T}) where T<:Union{Real,Complex}
    y = nfft(plan,reshape(x[:],shape))
    if symmetrize
      return vec(y) .*sqrt.(density)
    end
    return vec(y)
  end

  function ctprodu(y::Vector{T}) where T<:Union{Real,Complex}
    if symmetrize
      x = nfft_adjoint(plan, y[:] .*sqrt.(density))
    else
      x = nfft_adjoint(plan, y[:])
    end

    # if isCircular(tr)
    #   circularShutter!(x,1.0)
    # end

    return vec(x)
  end

  function invprodu(y::Vector{T}) where T<:Union{Real,Complex}
    if symmetrize
      x = nfft_adjoint(plan, y[:] .* sqrt.(density))
    else
      x = nfft_adjoint(plan, y[:] .* density) * sqrt(prod(shape))
    end
    return vec(x)
  end

  return NFFTOp{ComplexF64,Nothing,Function,Function}(size(nodes,2), prod(shape), false, false
            , produ
            , nothing
            , ctprodu
            , invprodu
            , density
            , isCircular(tr) )
end

function adjoint(op::NFFTOp{T}) where T
  return LinearOperator{T,Function,Nothing,Function}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod, nothing, op.prod)
end
