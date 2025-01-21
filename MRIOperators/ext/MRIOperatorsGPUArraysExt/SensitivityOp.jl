function MRIOperators.prod_smap!(y::vecT, smaps::matT, x::vecT, numVox, numChan, numContr=1) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}}
  x_ = reshape(x,numVox,numContr)
  y_ = reshape(y,numVox,numChan,numContr)

  @assert size(smaps) == (size(y_,1), size(y_,2))

  @kernel inbounds = true cpu = false function smap_kernel!(y, x, smaps)
    cart = @index(Global, Cartesian)
    y[cart] = x[cart[1], cart[3]] * smaps[cart[1], cart[2]]
  end
  kernel! = smap_kernel!(get_backend(y))
  kernel!(y, x, smaps; ndrange = size(y))
  return y
end

function MRIOperators.ctprod_smap!(y::vecT, smapsC::matT, x::vecT, numVox, numChan, numContr=1) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}}
  x_ = reshape(x,numVox,numChan,numContr)
  y_ = reshape(y,numVox,numContr)

  @assert size(smapsC) == (size(x_,1), size(x_,2))

  y_ .= 0
  # Inner loop to avoid race condition
  @kernel inbounds = true cpu = false function smap_kernel_adjoint!(y, x, smapsC, limit)
    cart = @index(Global, Cartesian)
    for i = 1:limit
      y[cart] += x[cart[1], i, cart[2]] * smapsC[cart[1],i]
    end
  end
  kernel! = smap_kernel_adjoint!(get_backend(y))
  kernel!(y, x, smapsC, numChan; ndrange = size(y))
  return y
end