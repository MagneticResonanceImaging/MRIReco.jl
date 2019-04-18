export loadImage, saveImage

function loadImage(filename::String)
  ni = niread(filename)
  pixelspacing = voxel_size(ni.header)
  return makeAxisArray(ni.raw, pixelspacing)
end

function saveImage(filename::String, im::AbstractArray{T,5}, makeAbs::Bool=false) where T
  voxel_size = ustrip.(uconvert.(u"m", pixelspacing(im)[1:3]))

  orientation = zeros(Float32, 3, 4)
  orientation[1,1] = orientation[2,2] = orientation[3,3] = 1.0

  ni = NIVolume(makeAbs ? abs.(im) : im; voxel_size=voxel_size,
                              orientation = orientation)
  # TODO add  dim_info=dim_info)
  niwrite(filename, ni)
end
