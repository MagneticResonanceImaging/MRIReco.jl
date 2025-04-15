function MRIOperators.produ!(out::vecTc, x::vecTc, shape::NTuple{2,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::matTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, matTc <: AbstractGPUArray{Tc}}

  factor = Tc(-2 * 1im * pi)
  limit = prod(shape)
  fill!(out, zero(Tc))
  @kernel inbounds = true cpu = false function produ_kernel!(out, x, shape, nodes, times, echoOffset, disturbanceTerm)
    k = @index(Global, Linear)
    for nx = 1:shape[1]
      for ny = 1:shape[2]
        phi = (nodes[1, k] * (nx - shape[1] / 2 - 1) + nodes[2, k] * (ny - shape[2] / 2 - 1))
        e = exp(factor * phi - (times[k] - echoOffset) * disturbanceTerm[nx, ny])
        out[k] += x[nx, ny] * e
      end
    end
  end
  kernel! = produ_kernel!(get_backend(out))
  kernel!(out, reshape(x, shape), shape, nodes, times, echoOffset, disturbanceTerm; ndrange = limit)

  return out
end

function MRIOperators.produ!(out::vecTc, x::vecTc, shape::NTuple{3,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::arrTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, arrTc <: AbstractGPUArray{Tc, 3}}

  factor = Tc(-2 * 1im * pi)
  limit = prod(shape)
  fill!(out, zero(Tc))
  @kernel inbounds = true cpu = false function produ_kernel!(out, x, shape, nodes, times, echoOffset, disturbanceTerm)
    k = @index(Global, Linear)
    for nx = 1:shape[1]
      for ny = 1:shape[2]
        for nz = 1:shape[3]
          phi = (nodes[1, k] * (nx - shape[1] / 2 - 1) + nodes[2, k] * (ny - shape[2] / 2 - 1) + nodes[3, k] * (nz - shape[3] / 2 - 1))
          e = exp(factor * phi - (times[k] - echoOffset) * disturbanceTerm[nx, ny, nz])
          out[k] += x[nx, ny, nz] * e
        end
      end
    end
  end
  kernel! = produ_kernel!(get_backend(out))
  kernel!(out, reshape(x, shape), shape, nodes, times, echoOffset, disturbanceTerm; ndrange = limit)

  return out
end

function MRIOperators.ctprodu!(out::vecTc, x::vecTc, shape::NTuple{2,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::matTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, matTc <: AbstractGPUArray{Tc}}

  factor = Tc(-2 * 1im * pi)
  limit = prod(shape)
  fill!(out, zero(Tc))
  @kernel inbounds = true cpu = false function ctprodu_kernel!(out, x, shape, nodes, times, echoOffset, disturbanceTerm)
    linearIndex = @index(Global, Linear)

    cartIndex = CartesianIndices(shape)[linearIndex]

    nx = cartIndex[1]
    ny = cartIndex[2]
    for k = 1:size(nodes, 2)
      phi = (nodes[1,k]*(nx-shape[1]/2-1)+nodes[2,k]*(ny-shape[2]/2-1))
      e = exp(factor * phi - (times[k]-echoOffset)*disturbanceTerm[cartIndex])
      out[cartIndex] += x[k] * conj(e)
    end
  end
  kernel! = ctprodu_kernel!(get_backend(out))
  kernel!(reshape(out, shape), x, shape, nodes, times, echoOffset, disturbanceTerm; ndrange = limit)

  return out
end

function MRIOperators.ctprodu!(out::vecTc, x::vecTc, shape::NTuple{3,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::arrTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, arrTc <: AbstractGPUArray{Tc, 3}}

  factor = Tc(-2 * 1im * pi)
  limit = prod(shape)
  fill!(out, zero(Tc))
  @kernel inbounds = true cpu = false function ctprodu_kernel!(out, x, shape, nodes, times, echoOffset, disturbanceTerm)
    linearIndex = @index(Global, Linear)

    cartIndex = CartesianIndices(shape)[linearIndex]

    nx = cartIndex[1]
    ny = cartIndex[2]
    nz = cartIndex[3]
    for k = 1:size(nodes, 2)
      phi = (nodes[1,k]*(nx-shape[1]/2-1)+nodes[2,k]*(ny-shape[2]/2-1)+nodes[3,k]*(nz-shape[3]/2-1))
      e = exp(factor * phi - (times[k]-echoOffset)*disturbanceTerm[cartIndex])
      out[cartIndex] += x[k] * conj(e)
    end
  end
  kernel! = ctprodu_kernel!(get_backend(out))
  kernel!(reshape(out, shape), x, shape, nodes, times, echoOffset, disturbanceTerm; ndrange = limit)
    
  return out
end