function MRIOperators.circularShutter!(I::AbstractGPUArray, radiusFactor::Number=1.0)
  imgSize = size(I)
  center = imgSize./2.0
  radius = maximum(center) * radiusFactor
  # Applying filtering
  gpu_call(I, center, radius) do ctx, I_, center_, radius_
    cart = @cartesianidx(I_)
    if sqrt(sum((Tuple(cart) .-center_).^2)) > radius_
        I[cart] = 0
    end
    return nothing
  end
  return I
end


function MRIOperators.circularShutterFreq!(I::AbstractGPUArray, radiusFactor::T=1.0) where T<:Number
  imgSize = size(I)
  fftI = fftshift(fft(I))
  center = imgSize./2.0
  radius = maximum(center) * radiusFactor
  # Applying filtering
  gpu_call(fftI, center, radius) do ctx, fftI_, center_, radius_
    cart = @cartesianidx(fftI_)
    if sqrt(sum((Tuple(cart) .-center_).^2)) > radius_
      fftI_[cart] = 0
    end
    return nothing
  end
  return MRIOperators.preserveType(I, ifft(ifftshift(fftI)))
end
