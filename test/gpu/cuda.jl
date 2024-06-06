using CUDA, CuNFFT

arrayTypes = [CuArray]

include(joinpath(@__DIR__(), "..", "runtests.jl"))