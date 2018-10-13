export circularShutterFreq,circularShutterSpat,getCircularMask,circularShutter

#TODO Remove when code has been cleaned up
function circularShutter(I::Matrix,radiusFactor::T=1.0) where T<:Number
  circularShutterFreq(I,radiusFactor)
end

function circularShutterFreq(I::Matrix,radiusFactor::T=1.0) where T<:Number
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
    return ifft(fftI)
end

function circularShutterFreq(I::Array{T,3},radiusFactor::T=1.0) where T<:Number
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
    return ifft(fftI)
end

function circularShutterSpat!(I::Matrix,radiusFactor::T=1.0) where T<:Number
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

function getCircularMask(shape::Tuple,radiusFactor::Float64=1.0)
  A = ones(Float64,shape)
  circularShutterSpat!(A,radiusFactor)
  return A
end
