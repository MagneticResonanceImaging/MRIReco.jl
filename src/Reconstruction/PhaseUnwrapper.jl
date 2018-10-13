
# Unwrap_phase


# http://scikit-image.org/docs/dev/auto_examples/plot_phase_unwrap.html
# Conda.add("scikit-image")
#=
using PyCall

@pyimport skimage.restoration as resto

function unwrap_phase(image::Matrix)
  resto.unwrap_phase(image)
end
=#
