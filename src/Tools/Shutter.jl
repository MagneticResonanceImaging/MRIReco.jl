export circularShutter!, circularShutterFreq!, getCircularMask

preserveType(orig::Array{T}, data::Array{U}) where {T<:Real, U<:Complex} = real(data)

preserveType(orig::Array{T}, data::Array{U}) where {T, U} = data

function circularShutter!(I::Matrix, radiusFactor::Number=1.0)
    imgSize = size(I)
    center = (imgSize[1]/2.0, imgSize[2]/2.0)
    radius = (center[1]>center[2] ? center[1] : center[2]) * radiusFactor
    # Applying filtering
    for j=1:imgSize[2]
        for i=1:imgSize[1]
            if sqrt((i-center[1])^2 + (j-center[2])^2) > radius
                I[i,j] = 0
            end
        end
    end
end

function getCircularMask(shape::Tuple, radiusFactor::Number=1.0)
  A = ones(Float64,shape)
  circularShutter!(A,radiusFactor)
  return A
end

function circularShutterFreq!(I::Matrix, radiusFactor::T=1.0) where T<:Number
    imgSize = size(I)
    fftI = fftshift(fft(I))
    center = (imgSize[1]/2.0, imgSize[2]/2.0)
    radius = (center[1]>center[2] ? center[1] : center[2]) * radiusFactor
    # Applying filtering
    for j=1:imgSize[2]
        for i=1:imgSize[1]
            if sqrt((i-center[1])^2 + (j-center[2])^2) > radius
                fftI[i,j] = 0
            end
        end
    end
    return preserveType(I, ifft(ifftshift(fftI)))
end

function circularShutterFreq!(I::Array{T,3}, radiusFactor::T=1.0) where T<:Number
    imgSize = size(I)
    fftI = fftshift(fft(I))
    center = (imgSize[1]/2.0, imgSize[2]/2.0, imgSize[3]/2.0)
    radius = (center[1]>center[2] ? center[1] : center[2]) * radiusFactor
    # Applying filtering
    for k=1:imgSize[3]
        for j=1:imgSize[2]
          for i=1:imgSize[1]
            if sqrt((i-center[1])^2 + (j-center[2])^2 + (k-center[3])^2) > radius
                fftI[i,j,k] = 0
            end
          end
        end
    end
    return preserveType(I, ifft(ifftshift(fftI)))
end
