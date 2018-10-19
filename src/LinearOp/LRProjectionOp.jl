export LRMapOp, LRProj!

mutable struct LRMapOp{T,F1<:FuncOrNothing,F2<:FuncOrNothing,F3<:FuncOrNothing} <:
                      AbstractLinearOperator{T,F1,F2,F3}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: F1
  tprod :: F2
  ctprod :: F3
  PCmap #:: Vector{Int}
end

function LRMapOp(phi, N)
  T,K,L = size(phi)

  nrow = K*L*N
  ncol = N*T

  phi_t =  permutedims( phi, [2,1,3] )

  PCmap = zeros(Int,N)

  function mul(x::Vector{T}) where T
    return LRMap(x,phi,PCmap)
  end
  function ctmul(y::Vector{T}) where T
    return invLRMap(y,phi_t,PCmap)
  end

  return LRMapOp{ComplexF64,Function,Nothing,Function}(K*L*N, N*T, false, false
                          , mul
                          , nothing
                          , ctmul
                          , PCmap )

end

function LRMap(x::Vector{V}, phi::Array{U,3}, PCmap::Vector{Int}) where {V,U}
  T,K,L = size(phi)                     # number of samples (echoes), number of PCS, number of subspaces
  N = floor(Int, length(x)/size(phi,1)) # number of signals (voxels)

  y = zeros(ComplexF64, N*K,L)

  PCmap[:] = segmentation( x, phi,N,size(phi,3) )

  for i = 1:L
    idx = findall(x->x==i, PCmap)
    fill!(y[:,i], 0)
    for j = 1:length(idx)
      y[idx[j]:N:end,i] = vec( transpose(x[idx[j]:N:end])*phi[:,:,i] )
    end
  end

  return vec(y)
end

function invLRMap(y::Vector{V}, phi_t::Array{U,3}, PCmap::Vector{Int}) where {V,U}

  K,T,L = size(phi_t)
  N = floor(Int, length(y)/(K*L))

  x = zeros(ComplexF64, N, T)
  for i = 1:L
    idx = find(x->x==i, PCmap)
    for j = 1:length(idx)
      x[idx[j]:N:end] += vec( transpose(y[(i-1)*N*K+idx[j]:N:i*N*K])*phi_t[:,:,i] )
    end
  end

  return vec(x)
end

function LRProj!(x, phi)
  T,K,L = size(phi)                     # number of samples (echoes), number of PCS, number of subspaces
  N = floor(Int, length(x)/size(phi,1)) # number of signals (voxels)
  phi² = zeros(eltype(phi),T,T,L)
  for i=1:L
    phi²[:,:,i] = phi[:,:,i]*phi[:,:,i]'
  end

  PCmap = segmentation( x, phi,N,size(phi,3) )
  for i=1:N
    x[i:N:end] = phi²[:,:,PCmap[i]]*x[i:N:end]
  end
end

#
# perform segmentation of image data x based on temporal subspaces phi
#
function segmentation(x::Vector{ComplexF64}, phi::Array{U,3}, N::Int64, L::Int64) where U
  dict = [phi[:,1,i] for i=1:L]
  subspaceMap = zeros(Int8, N)
  for i = 1:N
    subspaceMap[i] = matchSignal(x[i:N:end], dict)
  end
  return subspaceMap
end

#
# for a vector x find the element of dict, that has the largest overlap
#
function matchSignal(x::Vector{ComplexF64}, dict::Vector{Vector{U}}) where U
  overlap = zeros(length(dict))
  for i = 1:length(dict)
    overlap[i] = abs(dot(x, dict[i]))^2
  end
  coeff, idx = findmax(overlap)
  return idx
end
