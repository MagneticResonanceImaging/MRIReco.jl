"""
    makeAxisArray(I::AbstractArray{T,5}, acqData::AcquisitionData) where T

creates an axies array with properly characterized axis from the image `I`.
For this, it uses the information in `acqData`.
The axies respectively describe the following coordinates:
x, y, z, acqData.numEchoes, acqData.numCoils
"""
function ImageUtils.makeAxisArray(I::AbstractArray{T,5}, acqData::AcquisitionData) where T

  offset = [0.0, 0.0, 0.0]*Unitful.mm
  spacing = pixelspacing(acqData)

  im = AxisArray(I,
		   Axis{:x}(range(offset[1], step=spacing[1], length=size(I,1))),
		   Axis{:y}(range(offset[2], step=spacing[2], length=size(I,2))),
		   Axis{:z}(range(offset[3], step=spacing[3], length=size(I,3))),
		   Axis{:echos}(1:size(I,4)),
       Axis{:coils}(1:size(I,5)))

  #imMeta = ImageMeta(im, Dict{String,Any}())
  #return imMeta
  return im
end

function ImageUtils.makeAxisArray(I::AbstractArray{T,5}, spacing::Vector{Float64}) where T

  offset = [0.0, 0.0, 0.0]*Unitful.mm

  sp = uconvert.(Unitful.mm, spacing*Unitful.m)

  im = AxisArray(I,
		   Axis{:x}(range(offset[1], step=sp[1], length=size(I,1))),
		   Axis{:y}(range(offset[2], step=sp[2], length=size(I,2))),
		   Axis{:z}(range(offset[3], step=sp[3], length=size(I,3))),
		   Axis{:echos}(1:size(I,4)),
       Axis{:coils}(1:size(I,5)))

  return im
end
