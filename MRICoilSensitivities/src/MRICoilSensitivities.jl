module MRICoilSensitivities

export estimateCoilSensitivities, mergeChannels, espirit, geometricCC_2d

using MRIBase
using FLoops
using LinearAlgebra
using FFTW

include("simple.jl")
include("espirit.jl")
include("CoilCompression.jl")

end # module
