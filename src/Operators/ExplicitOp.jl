export ExplicitOp

mutable struct ExplicitOp{T,F1<:FuncOrNothing,F2<:FuncOrNothing,F3<:FuncOrNothing} <: AbstractLinearOperator{T,F1,F2,F3}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: F1
  ctprod :: F2
  inv :: F3
  density::Vector{Float64}
end

#
# Linear Operator to perform explicite reconstruction
#
function ExplicitOp(shape::NTuple{D,Int64}, tr::Trajectory, correctionmap::Array{ComplexF64,D}
                        ; symmetrize::Bool=true
                        , echoImage::Bool=false) where D

  nodes,times = kspaceNodes(tr), readoutTimes(tr)
  nrow = size(nodes,2)
  ncol = prod(shape)

  density = zeros(size(nodes,2))
  if symmetrize
      plan = NFFTPlan(nodes, shape, 4, 1.75)
      density = convert(Vector{Float64}, sdc(plan))
  end

  # if echo image is desired echo time is needed as an offset
  if echoImage
      echoOffset = echoTime(tr)
  else
      echoOffset = 0.0
  end

  return ExplicitOp{ComplexF64,Nothing,Function,Function}(nrow, ncol, false, false
            , x->produ(x, shape, nodes, times, echoOffset, correctionmap, density, symmetrize)
            , nothing
            , y->ctprodu(y, shape, nodes, times, echoOffset, correctionmap, density, symmetrize)
            , y->ctprodu(y, shape, nodes, times, echoOffset, correctionmap, density, symmetrize)
            , density )
end

function produ(x::Vector{T}, shape::NTuple{2,Int64},
                    nodes::Matrix{Float64}, times::Vector{Float64}, echoOffset::Float64,
                    correctionmap::Matrix{ComplexF64}, density::Vector{Float64},
                    symmetrize::Bool) where T<:Union{Real,Complex}
   if symmetrize
       x = x .* sqrt(density)
   end

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
                    correctionmap::Array{ComplexF64,3}, density::Vector{Float64},
                    symmetrize::Bool) where T<:Union{Real,Complex}
   if symmetrize
       x = x .* sqrt(density)
   end

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
                        correctionmap::Matrix{ComplexF64}, density::Vector{Float64},
                        symmetrize::Bool) where T<:Union{Real,Complex}

  if symmetrize
      x = x .* sqrt(density) # <- using for FISTA
  end

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
                        correctionmap::Array{ComplexF64,3}, density::Vector{Float64},
                        symmetrize::Bool) where T<:Union{Real,Complex}

  if symmetrize
      x = x .* sqrt(density) # <- using for FISTA
  end

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
  return LinearOperator{T,Function,Nothing,Function}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod, nothing, op.prod)
end
