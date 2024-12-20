using CairoMakie
function exportImage(filepath::AbstractString,img)
  f = Figure(figure_padding = 0)
  ax=Axis(f[1,1],aspect=1)
  heatmap!(ax,rotr90(img ),colormap=:grays)
  hidedecorations!(ax)
  colsize!(f.layout, 1, Aspect(1,1))
  resize_to_layout!(f)
  f
  save(filepath,f)
end