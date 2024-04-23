function MRIOperators.prod_smap!(y::vecT, smaps::matT, x::vecT, numVox, numChan, numContr=1) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}}
  x_ = reshape(x,numVox,numContr)
  y_ = reshape(y,numVox,numChan,numContr)

  @assert size(smaps) == (size(y_,1), size(y_,2))

  gpu_call(y_, x_, smaps) do ctx, ygpu, xgpu, smapsgpu
    cart = @cartesianidx(ygpu)
    ygpu[cart] = xgpu[cart[1],cart[3]] * smapsgpu[cart[1],cart[2]]
    return nothing
  end
  return y
end

function MRIOperators.ctprod_smap!(y::vecT, smapsC::matT, x::vecT, numVox, numChan, numContr=1) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}}
  x_ = reshape(x,numVox,numChan,numContr)
  y_ = reshape(y,numVox,numContr)

  @assert size(smapsC) == (size(x_,1), size(x_,2))

  y_ .= 0
  # Inner loop to avoid race condition
  gpu_call(y_, x_, smapsC, numChan) do ctx, ygpu, xgpu, smapsCgpu, limit
    cart = @cartesianidx(ygpu)
    for i = 1:limit
      ygpu[cart] += xgpu[cart[1], i, cart[2]] * smapsCgpu[cart[1],i]
    end
    return nothing
  end
  return y
end