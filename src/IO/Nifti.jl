export loadImage, saveImage, makeImage

function loadImage(filename::String)
  ni = niread(filename)
  pixelspacing = voxel_size(ni.header)
  return makeAxisArray(ni.raw, pixelspacing)
end

function saveImage(filename::String, im::AbstractArray{T,5}, makeAbs::Bool=false) where T
  voxel_size = ustrip.(uconvert.(u"m", pixelspacing(im)[1:3]))

  ni = NIVolume(makeAbs ? abs.(im) : im; voxel_size=voxel_size)
  # TODO add orientation   orientation=orientation, dim_info=dim_info)
  niwrite(filename, ni)
end
