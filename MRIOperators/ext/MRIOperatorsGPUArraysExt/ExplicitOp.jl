function MRIOperators.produ!(out::vecTc, x::vecTc, shape::NTuple{2,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::matTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, matTc <: AbstractGPUArray{Tc}}

  factor = Tc(-2 * 1im * pi)
  limit = prod(shape)
  gpu_call(out, reshape(x, shape), shape, nodes, times, echoOffset, disturbanceTerm; elements = size(nodes, 2)) do ctx, out_, x_, shape_, nodes_, times_, echoOffset_, disturbanceTerm_
    k = linear_index(ctx)
    k = linear_index(ctx)
    if !(1 <= k <= limit)
      return nothing
    end

    for nx = 1:shape_[1]
      for ny = 1:shape_[2]
        phi = (nodes_[1, k] * (nx - shape_[1] / 2 - 1) + nodes_[2, k] * (ny - shape_[2] / 2 - 1))
        e = exp(factor * phi - (times_[k] - echoOffset_) * disturbanceTerm_[nx, ny])
        out_[k] += x_[nx, ny] * e
      end
    end
    return nothing
  end

  return out
end

function MRIOperators.produ!(out::vecTc, x::vecTc, shape::NTuple{3,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::arrTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, arrTc <: AbstractGPUArray{Tc, 3}}

  factor = Tc(-2 * 1im * pi)
  limit = prod(shape)
  gpu_call(out, reshape(x, shape), shape, nodes, times, echoOffset, disturbanceTerm; elements = size(nodes, 2)) do ctx, out_, x_, shape_, nodes_, times_, echoOffset_, disturbanceTerm_
    k = linear_index(ctx)
    if !(1 <= k <= limit)
      return nothing
    end

    for nx = 1:shape_[1]
      for ny = 1:shape_[2]
        for nz = 1:shape_[3]
          phi = (nodes_[1, k] * (nx - shape_[1] / 2 - 1) + nodes_[2, k] * (ny - shape_[2] / 2 - 1) + nodes_[3, k] * (nz - shape_[3] / 2 - 1))
          e = exp(factor * phi - (times_[k] - echoOffset_) * disturbanceTerm_[nx, ny, nz])
          out_[k] += x_[nx, ny, nz] * e
        end
      end
    end
    return nothing
  end

  return out
end

function MRIOperators.ctprodu!(out::vecTc, x::vecTc, shape::NTuple{2,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::matTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, matTc <: AbstractGPUArray{Tc}}

  out = reshape(out, shape)
  for nx = 1:shape[1]
    for ny = 1:shape[2]
      for k = 1:size(nodes, 2)
        phi = (nodes[1, k] * (nx - shape[1] / 2 - 1) + nodes[2, k] * (ny - shape[2] / 2 - 1))
        e = exp(-2 * 1im * pi * phi - (times[k] - echoOffset) * disturbanceTerm[nx, ny])
        out[nx, ny] += x[k] * conj(e)
      end
    end
  end
  return out
end

function MRIOperators.ctprodu!(out::vecTc, x::vecTc, shape::NTuple{3,Int64},
  nodes::matT, times::vecT, echoOffset::T,
  disturbanceTerm::arrTc) where {T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractGPUArray{Tc}, matT <: AbstractGPUArray{T}, vecT <: AbstractGPUArray{T}, arrTc <: AbstractGPUArray{Tc, 3}}

  out = reshape(out, shape)
  for nx = 1:shape[1]
    for ny = 1:shape[2]
      for nz = 1:shape[3]
        for k = 1:size(nodes, 2)
          phi = (nodes[1, k] * (nx - shape[1] / 2 - 1) + nodes[2, k] * (ny - shape[2] / 2 - 1) + nodes[3, k] * (nz - shape[3] / 2 - 1))
          e = exp(-2 * 1im * pi * phi - (times[k] - echoOffset) * disturbanceTerm[nx, ny, nz])
          out[nx, ny, nz] += x[k] * conj(e)
        end
      end
    end
  end
  return out
end