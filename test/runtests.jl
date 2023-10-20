using MRIReco, MRISampling, MRISimulation, FFTW
using Test
using LinearAlgebra
using Random
using ImageUtils
using Scratch

const tmpdir  = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpdir."

include("testIO.jl")
include("testReconstruction.jl")
include("testSpecificApplications.jl")
