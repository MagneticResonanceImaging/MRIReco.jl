export RecoFileIBI, recoImage, recoParams, saveasRecoFile

struct RecoFileIBI <: MRIFile
  filename::String
end

function recoImage(f::RecoFileIBI)
  img = h5read(f.filename, "/image")
  # workaround for hdf5 not supporting complex
  imgComplex = collect(reinterpret(Complex{eltype(img)}, vec(img)))

  return reshape(imgComplex, (size(img)[2:end]...,))
end

function recoParams(f::RecoFileIBI)
  return loadParams(f.filename, "params")
end

function saveasRecoFile(filename::AbstractString, img::Array{Complex{T}}, params::Dict) where T<:Real
  h5open(filename, "w") do file
    # img_real = reshape(reinterpret(T, vec(img)), (2,size(img)...,)) # workaround for hdf5 not supporting complex
    img_real = reinterpret(T, vec(img))
    img_real = collect(reshape(img_real,(2,size(img)...,)))
    write(file, "/image", img_real)
    saveParams(file, "params", params)
  end
end
