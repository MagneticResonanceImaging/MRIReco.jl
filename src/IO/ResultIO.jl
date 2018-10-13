export saveImageData, loadImageData

#FIXME: Do we really need this? => check if nifti does not cover it.

function saveImagesResults(recoImages::Vector{Vector{Float64}},shape::Tuple{Int64,Int64},imageSavePath::String,skip::Vector{Int64})
    count = 1
    for i=1:length(recoImages)
      while skip[count] == 1
          count+=1
      end
      saveString = string(imageSavePath,"$count/corrected.png")
      im = abs(reshape(recoImages[i],shape))
      im = colorview(Gray,im./maximum(im))
      save(saveString,im)
      count+=1
    end
end

function saveImagesResults(recoImages::Vector{Vector{Float64}},shape::Tuple{Int64,Int64,Int64},imageSavePath::String,skip::Vector{Int64})
    count = 1
    for i=1:length(recoImages)
      while skip[count] == 1
          count+=1
      end
      saveString = string(imageSavePath,"$count/corrected.nrrd")
      save(saveString,reshape(recoImages[i],shape))
      saveString = string(imageSavePath,"$count/correctedSlice.png")
      im = abs(reshape(recoImages[i],shape))
      im = colorview(Gray,im[round(Int64,shape[1]/2),:,:]./maximum(im[round(Int64,shape[1]/2),:,:]))
      save(saveString,im)
      count+=1
    end
end


#
# save image data. Images are stored in 3d-Arrays:
# 1. dimension: spatial
# 2. dimension: temporal
# 3. dimension: receive coils
#
function saveImageData(filename::AbstractString, image::Array{Float64,3}, shape::Tuple)

  h5open(filename, "w") do file

    write(file, "shape", collect(shape)) # shape for the spatial dimensions
    write(file, "image", image)
  end
end

#
# returns both the 3d Array containing image data and the spatial shape of the image
#
function loadImageData(filename::AbstractString)
  shape = Tuple( h5read(filename, "shape") )
  image = h5read(filename, "image")

  return image, shape
end
