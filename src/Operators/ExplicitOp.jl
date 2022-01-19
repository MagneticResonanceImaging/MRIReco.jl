export ExplicitOp

mutable struct ExplicitOp{T,F1,F2} <: AbstractLinearOperator{T}
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
  Mv5 :: Vector{T}
  Mtu5 :: Vector{T}
end

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
function ExplicitOp(shape::NTuple{D,Int64}, tr::Trajectory, correctionmap::Array{ComplexF64,D}
                        ; echoImage::Bool=false, kargs...) where D

  nodes,times = kspaceNodes(tr), readoutTimes(tr)
  nrow = size(nodes,2)
  ncol = prod(shape)

  # if echo image is desired echo time is needed as an offset
  if echoImage
      echoOffset = echoTime(tr)
  else
      echoOffset = 0.0
  end

  return ExplicitOp{ComplexF64,Nothing,Function}(nrow, ncol, false, false
            , (res,x)->(res .= produ(x, shape, nodes, times, echoOffset, correctionmap))
            , nothing
            , (res,y)->(res .= ctprodu(y, shape, nodes, times, echoOffset, correctionmap))
            , 0,0,0, false, false, false, ComplexF64[], ComplexF64[])
end

function produ(x::Vector{T}, shape::NTuple{2,Int64},
                    nodes::Matrix{Float64}, times::Vector{Float64}, echoOffset::Float64,
                    correctionmap::Matrix{ComplexF64}) where T<:Union{Real,Complex}

   if isempty(correctionmap)
       disturbanceTerm = zeros(ComplexF64,shape...)
   else
       disturbanceTerm = correctionmap
   end

   x= reshape(x,shape)
   out = zeros(ComplexF64,size(nodes,2))
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
   return vec(out)
end

function produ(x::Vector{T}, shape::NTuple{3,Int64},
                    nodes::Matrix{Float64}, times::Vector{Float64}, echoOffset::Float64,
                    correctionmap::Array{ComplexF64,3}) where T<:Union{Real,Complex}

   if isempty(correctionmap)
       disturbanceTerm = zeros(ComplexF64,shape...)
   else
       disturbanceTerm = correctionmap
   end

   x = reshape(x,shape)
   out = zeros(ComplexF64,size(nodes,2))
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
   return vec(out)
end

function ctprodu(x::Vector{T}, shape::NTuple{2,Int64},
                        nodes::Matrix{Float64}, times::Vector{Float64}, echoOffset::Float64,
                        correctionmap::Matrix{ComplexF64}) where T<:Union{Real,Complex}

  if isempty(correctionmap)
      disturbanceTerm = zeros(ComplexF64,shape...)
  else
      disturbanceTerm = correctionmap
  end

  out = zeros(ComplexF64,shape)
  for nx=1:shape[1]
      for ny=1:shape[2]
          for k=1:size(nodes,2)
              phi = (nodes[1,k]*(nx-shape[1]/2-1)+nodes[2,k]*(ny-shape[2]/2-1))
              e = exp(-2*1im*pi*phi - (times[k]-echoOffset)*disturbanceTerm[nx,ny])
              out[nx,ny] += x[k] * conj(e)
          end
      end
  end
  return vec(out)
end

function ctprodu(x::Vector{T}, shape::NTuple{3,Int64},
                        nodes::Matrix{Float64}, times::Vector{Float64}, echoOffset::Float64,
                        correctionmap::Array{ComplexF64,3}) where T<:Union{Real,Complex}

  if isempty(correctionmap)
      disturbanceTerm = zeros(ComplexF64,shape...)
  else
      disturbanceTerm = correctionmap
  end

  out = zeros(ComplexF64,shape)
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
  return vec(out)
end

function adjoint(op::ExplicitOp{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!)
end
