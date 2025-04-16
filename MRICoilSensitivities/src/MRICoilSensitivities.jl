module MRICoilSensitivities

export estimateCoilSensitivities, mergeChannels, espirit, geometricCoilCompression

using MRIBase
using OhMyThreads: @tasks, @local
using LinearAlgebra
using FFTW

include("simple.jl")
include("espirit.jl")
include("CoilCompression.jl")

end # module
