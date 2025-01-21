@kernel inbounds = true cpu = false function circularShutterKernel(I, center, radius)
  cart = @index(Global, Cartesian)
  if sqrt(sum((Tuple(cart) .-center).^2)) > radius
    I[cart] = 0
  end
end

function MRIOperators.circularShutter!(I::AbstractGPUArray, radiusFactor::Number=1.0)
  imgSize = size(I)
  center = imgSize./2.0
  radius = maximum(center) * radiusFactor
  kernel!= circularShutterKernel(get_backend(I))
  kernel!(I, center, radius; ndrange = imgSize)
  return I
end



function MRIOperators.circularShutterFreq!(I::AbstractGPUArray, radiusFactor::T=1.0) where T<:Number
  imgSize = size(I)
  fftI = fftshift(fft(I))
  center = imgSize./2.0
  radius = maximum(center) * radiusFactor
  kernel!= circularShutterKernel(get_backend(I))
  kernel!(fftI, center, radius; ndrange = size(fftI))
  return MRIOperators.preserveType(I, ifft(ifftshift(fftI)))
end
