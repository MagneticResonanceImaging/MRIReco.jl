export translatePhantom

function translatePhantom(A::Matrix{T},x::Float64=0.,y::Float64=0.) where T<:Real

freq = fftshift(fft(A))

posY = collect((0:size(freq,1)-1)/size(freq,1) .- 0.5)
posX = collect((0:size(freq,2)-1)/size(freq,2) .- 0.5)

for xInd=1:size(freq,2)
  for yInd=1:size(freq,1)
    freq[yInd,xInd] *= exp(-1im*2*pi*(posY[yInd]*y+posX[xInd]*x))
  end
end

return real(ifft(ifftshift(freq)))

end

function translatePhantom(A::Array{T,3},x::Float64=0.,y::Float64=0.,z::Float64=0.) where T<:Real

freq = fftshift(fft(A))

posY = collect((0:size(freq,1)-1)/size(freq,1) .- 0.5)
posX = collect((0:size(freq,2)-1)/size(freq,2) .- 0.5)
posZ = collect((0:size(freq,3)-1)/size(freq,3) .- 0.5)

for zInd=1:size(freq,3)
  for xInd=1:size(freq,2)
    for yInd=1:size(freq,1)
      freq[yInd,xInd,zInd] *= exp(-1im*2*pi*(posY[yInd]*y+posX[xInd]*x+posZ[zInd]*z))
    end
  end
end

return real(ifft(ifftshift(freq)))

end
