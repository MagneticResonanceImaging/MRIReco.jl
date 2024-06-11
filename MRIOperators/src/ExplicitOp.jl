export ExplicitOp

mutable struct ExplicitOp{T, vecT <: AbstractVector{T}, F1, F2} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod! :: Function
  tprod! :: F1
  ctprod! :: F2
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  args5 :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5 :: vecT
  Mtu5 :: vecT
end

LinearOperators.storage_type(op::ExplicitOp) = typeof(op.Mv5)

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
"""
function ExplicitOp(shape::NTuple{D,Int64}, tr::Trajectory{T}, correctionmap::AbstractArray{Tc,D}
                        ; echoImage::Bool=false, S = storage_type(correctionmap), kargs...) where {T, Tc <: Union{Complex{T}, T}, D}

  nodes,times = kspaceNodes(tr), readoutTimes(tr)
  nrow = size(nodes,2)
  ncol = prod(shape)

  if isempty(correctionmap)
    disturbanceTerm = zeros(Complex{T}, shape...)
  else
    disturbanceTerm = correctionmap
  end

  tmp = S(undef, 0)
  if !isa(tmp, Vector)
    nodes = copyto!(similar(tmp, T, size(nodes)), nodes)
    times = copyto!(similar(tmp, T, size(times)), times)
    disturbanceTerm = copyto!(similar(tmp, size(disturbanceTerm)), disturbanceTerm)
  end

  # if echo image is desired echo time is needed as an offset
  if echoImage
      echoOffset = echoTime(tr)
  else
      echoOffset = 0.0
  end
  echoOffset = T(echoOffset)

  return ExplicitOp{Complex{T}, S, Nothing, Function}(nrow, ncol, false, false
            , (res,x)->(produ!(res, x, shape, nodes, times, echoOffset, disturbanceTerm))
            , nothing
            , (res,y)->(ctprodu!(res, y, shape, nodes, times, echoOffset, disturbanceTerm))
            , 0,0,0, false, false, false, tmp, tmp)
end

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

function Base.adjoint(op::ExplicitOp{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!; S = LinearOperators.storage_type(op))
end
