"""
    makeAxisArray(I::AbstractArray{T,5}, acqData::AcquisitionData) where T

creates an axes array with properly characterized axis from the image `I`.
For this, it uses the information in `acqData`.
The axies respectively describe the following coordinates:
x, y, z, acqData.numEchoes, acqData.numCoils
"""
function makeAxisArray(I::AbstractArray{T,6}, acqData::AcquisitionData) where T

  offset = [0.0, 0.0, 0.0]*Unitful.mm
  spacing = fieldOfView(acqData)./encodingSize(acqData)*Unitful.mm

  im = AxisArray(I,
		   Axis{:x}(range(offset[1], step=spacing[1], length=size(I,1))),
		   Axis{:y}(range(offset[2], step=spacing[2], length=size(I,2))),
		   Axis{:z}(range(offset[3], step=spacing[3], length=size(I,3))),
		   Axis{:echos}(1:size(I,4)),
       Axis{:coils}(1:size(I,5)),
       Axis{:repetitions}(1:size(I,6)))

  #imMeta = ImageMeta(im, Dict{String,Any}())
  #return imMeta
  return im
end
