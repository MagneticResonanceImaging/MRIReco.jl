using FileIO

export exportImage

function exportImage(filename, data_::Matrix)
 file, ext = splitext(filename)
 filename_ = file*".png"

 #data  = squeeze(data_)

  data_ = repeat(data_,inner=[1,1])

 argbdata = colorize( data_, 0, 1, cmap("gray"), normalize=true)
 rgbdata = convert(Array{RGB},argbdata)
 #rgbdata = flipdim(rgbdata,1)

 save(filename_, rgbdata)
end


function exportImage(filename, data_::Matrix, vmin, vmax,colormap="gray")
 file, ext = splitext(filename)
 filename_ = file*".png"

 #data  = squeeze(data_)

  data_ = repeat(data_,inner=[1,1])

 argbdata = colorize( data_, vmin, vmax, cmap(colormap), normalize=false )
 rgbdata = convert(Array{RGB},argbdata)
 #rgbdata = flipdim(rgbdata,1)

 save(filename_, rgbdata)
end




