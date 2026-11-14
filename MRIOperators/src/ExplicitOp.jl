export ExplicitOp

mutable struct ExplicitOp{T, vecT <: AbstractVector{T}, F1, F2, D, N, Ti, R, M} <: AbstractMRIOperator{T}
  const nrow :: Int
  const ncol :: Int
  const symmetric :: Bool
  const hermitian :: Bool
  const prod! :: Function
  const tprod! :: F1
  const ctprod! :: F2
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  Mv :: vecT
  Mtu :: vecT
  # inputs of prod! and ctprod!, kept for Adapt
  const shape :: NTuple{D,Int}
  const nodes :: N
  const times :: Ti
  const echoOffset :: R
  const disturbanceTerm :: M
end

LinearOperators.storage_type(::ExplicitOp{T, vecT}) where {T, vecT} = vecT

"""
    ExplicitOp(shape::NTuple{D,Int64}, tr::Trajectory, correctionmap::Array{ComplexF64,D}
            ; echoImage::Bool=false, kargs...) where D

generates a `ExplicitOp` which explicitely evaluates the MRI Fourier signal encoding operator.

# Arguments:
* `shape::NTuple{D,Int64}`             - size of image to encode/reconstruct
* `tr::Trajectory`                     - Trajectory with the kspace nodes to sample
* `correctionmap::Array{ComplexF64,D}` - fieldmap for the correction of off-resonance effects
* `echoImage::Bool=false`              - if true sampling times will only be considered relative to the echo time
                                         this results in complex valued image even for real-valued input.
* `kargs`                              - additional keyword arguments

# Moving to another array backend

    Adapt.adapt(CuArray, op::ExplicitOp)

returns a new `ExplicitOp` built from `Adapt.adapt(CuArray, ...)` of the k-space nodes, readout times and
correction map of `op`, with the same `shape` and echo time offset; the trajectory is not used again.
`adjoint(op)` is an `AdjointLinearOperator` with parent `op`, and adapting it adapts that parent.
On the same backend, the arrays may be shared with `op` and are only read. `op` itself is not modified.
`SensitivityOp`, `SubspaceOp`, `MapSliceOp`, `ExplicitOp` and `FieldmapNFFTOp` are the MRI operators that support `Adapt.adapt` so far.
"""
function ExplicitOp(shape::NTuple{D,Int64}, tr::Trajectory{T}, correctionmap::AbstractArray{Tc,D}
                        ; echoImage::Bool=false, S = Vector{Complex{T}}, kargs...) where {T, Tc <: Union{Complex{T}, T}, D}

  nodes,times = kspaceNodes(tr), readoutTimes(tr)

  if isempty(correctionmap)
    disturbanceTerm = zeros(Complex{T}, shape...)
  else
    disturbanceTerm = correctionmap
  end

  baseArrayType = stripParameters(S)
  nodes = adapt(baseArrayType, nodes)
  times = adapt(baseArrayType, times)
  disturbanceTerm = adapt(baseArrayType, disturbanceTerm)

  # if echo image is desired echo time is needed as an offset
  if echoImage
      echoOffset = echoTime(tr)
  else
      echoOffset = 0.0
  end
  echoOffset = T(echoOffset)

  return ExplicitOp(Complex{T}, S, shape, nodes, times, echoOffset, disturbanceTerm)
end

# Builds the operator from the inputs of prod! and ctprod!; used by the constructor above and Adapt.
function ExplicitOp(::Type{Tc}, S, shape, nodes, times, echoOffset, disturbanceTerm) where Tc
  tmp = S(undef, 0)
  return ExplicitOp{Tc, S, Nothing, Function, length(shape), typeof(nodes), typeof(times), typeof(echoOffset), typeof(disturbanceTerm)}(
              size(nodes,2), prod(shape), false, false
            , (res,x)->(produ!(res, x, shape, nodes, times, echoOffset, disturbanceTerm))
            , nothing
            , (res,y)->(ctprodu!(res, y, shape, nodes, times, echoOffset, disturbanceTerm))
            , 0,0,0, tmp, tmp, shape, nodes, times, echoOffset, disturbanceTerm)
end

Adapt.adapt_structure(to, op::ExplicitOp{T}) where T =
  ExplicitOp(T, typeof(adapt(to, LinearOperators.storage_type(op)(undef, 0))), op.shape, adapt(to, op.nodes),
             adapt(to, op.times), op.echoOffset, adapt(to, op.disturbanceTerm))

# The adjoint of an MRI operator is adapted through its parent.
Adapt.adapt_structure(to, op::LinearOperators.AdjointLinearOperator{T, <:AbstractMRIOperator}) where T =
  adjoint(adapt(to, op.parent))

function produ!(out::Vector{Tc}, x::Vector{Tc}, shape::NTuple{2,Int64},
                    nodes::Matrix{T}, times::Vector{T}, echoOffset::T,
                    disturbanceTerm::Matrix{Tc}) where {T, Tc <: Union{T, Complex{T}}}

   x= reshape(x,shape)
   fill!(out, zero(Tc))
   for nx=1:shape[1]
       for ny=1:shape[2]
           for k=1:size(nodes,2)
               phi = (nodes[1,k]*(nx-shape[1]/2-1)+
                      nodes[2,k]*(ny-shape[2]/2-1))
               e = exp(-2*1im*pi*phi - (times[k]-echoOffset)*disturbanceTerm[nx,ny])
               out[k] += x[nx,ny] * e
           end
       end
   end
   return out
end

function produ!(out::Vector{Tc}, x::Vector{Tc}, shape::NTuple{3,Int64},
                    nodes::Matrix{T}, times::Vector{T}, echoOffset::T,
                    disturbanceTerm::Array{Tc,3}) where {T, Tc <: Union{T, Complex{T}}}


   x = reshape(x,shape)
   fill!(out, zero(Tc))
   for nx=1:shape[1]
       for ny=1:shape[2]
           for nz=1:shape[3]
               for k=1:size(nodes,2)
                   phi = (nodes[1,k]*(nx-shape[1]/2-1)+nodes[2,k]*(ny-shape[2]/2-1)+nodes[3,k]*(nz-shape[3]/2-1))
                   e = exp(-2*1im*pi*phi - (times[k]-echoOffset)*disturbanceTerm[nx,ny,nz])
                   out[k] += x[nx,ny,nz] * e
               end
           end
       end
   end
   return out
end

function ctprodu!(out::Vector{Tc}, x::Vector{Tc}, shape::NTuple{2,Int64},
                        nodes::Matrix{T}, times::Vector{T}, echoOffset::T,
                        disturbanceTerm::Matrix{Tc}) where {T, Tc <: Union{T, Complex{T}}}


  out = reshape(out, shape)
  fill!(out, zero(Tc))
  for nx=1:shape[1]
      for ny=1:shape[2]
          for k=1:size(nodes,2)
              phi = (nodes[1,k]*(nx-shape[1]/2-1)+nodes[2,k]*(ny-shape[2]/2-1))
              e = exp(-2*1im*pi*phi - (times[k]-echoOffset)*disturbanceTerm[nx,ny])
              out[nx,ny] += x[k] * conj(e)
          end
      end
  end
  return out
end

function ctprodu!(out::Vector{Tc}, x::Vector{Tc}, shape::NTuple{3,Int64},
                        nodes::Matrix{T}, times::Vector{T}, echoOffset::T,
                        disturbanceTerm::Array{Tc,3}) where {T, Tc <: Union{T, Complex{T}}}

  out = reshape(out, shape)
  fill!(out, zero(Tc))
  for nx=1:shape[1]
      for ny=1:shape[2]
          for nz=1:shape[3]
              for k=1:size(nodes,2)
                 phi = (nodes[1,k]*(nx-shape[1]/2-1)+nodes[2,k]*(ny-shape[2]/2-1)+nodes[3,k]*(nz-shape[3]/2-1))
                 e = exp(-2*1im*pi*phi - (times[k]-echoOffset)*disturbanceTerm[nx,ny,nz])
                 out[nx,ny,nz] += x[k] * conj(e)
              end
          end
      end
  end
  return out
end
