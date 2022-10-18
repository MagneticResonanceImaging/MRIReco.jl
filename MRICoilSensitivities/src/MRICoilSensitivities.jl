module MRICoilSensitivities

export estimateCoilSensitivities, mergeChannels, espirit, geometricCoilCompression

using MRIBase
using FLoops
using LinearAlgebra
using FFTW

include("simple.jl")
include("espirit.jl")
include("CoilCompression.jl")

end # module
