export relaxationMap

"""
  Calculates the relaxationMap based on the pixel values of the image

  ...
  # Arguments
  * `image::Matrix` : Image whose values determine the relaxationrate at each pixel.
  * `minRelaxation::Float64=5.0` : minimum relaxationrate measured in [Hz]
  * `maxRelaxation::Float64=125.0` : maximum relaxationrate measured in [Hz]
  ...
"""
function relaxationMap(image::Matrix,minRelaxation::Float64=5.0,maxRelaxation::Float64=125.0)
  # consistency check
  if minRelaxation > maxRelaxation
    warning("Minimal frequency > Maximalfrequency! Both will be swapped now")
    minRelaxation,maxRelaxation = maxRelaxation,minRelaxation
  end

  relaxMap = zeros(size(image))
  # normalize first
  minVal = minimum(image)
  maxVal = maximum(image)
  z = image .- minVal ./ (maxVal-minVal)

  # now denormalize
  relaxMap = z .* (maxRelaxation - minRelaxation) .+ minRelaxation

  return relaxMap

end
