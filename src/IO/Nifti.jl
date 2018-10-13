using NIfTI, Images

export loadImage, saveImage, makeImage

function loadImage(filename::String)
  ni = niread(filename)
  pixelspacing = voxel_size(ni.header)
  return  ni.raw #makeImage(ni.raw, pixelspacing)
end

#=
function makeImage{T}(a::AbstractArray{T,3}, pixelspacing)
  dims = (:x,:y,:z, :coils, :slices)
  im = AxisArray(a, (:x,:y,:z),
                      tuple(pixelspacing...),
                      tuple(0.0, 0.0, 0.0))

  imMeta = ImageMeta(im, Dict{String,Any}())
  return imMeta
end
=#

function saveImage(filename::String, im::AbstractArray)
  voxel_size = pixelspacing(im)

  ni = NIVolume(im)#; voxel_size=voxel_size)
  # TODO add orientation   orientation=orientation, dim_info=dim_info)
  niwrite(filename, ni)
end
