export loadImage, saveImage, makeImage

function loadImage(filename::String)
  ni = niread(filename)
  pixelspacing = voxel_size(ni.header)
  return  ni.raw #makeImage(ni.raw, pixelspacing)
end

function saveImage(filename::String, im::AbstractArray)
  voxel_size = pixelspacing(im)

  ni = NIVolume(im)#; voxel_size=voxel_size)
  # TODO add orientation   orientation=orientation, dim_info=dim_info)
  niwrite(filename, ni)
end
