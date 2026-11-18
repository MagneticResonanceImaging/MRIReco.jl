using Metal
using MRIReco, MRISampling, MRISimulation, FFTW
using MRIReco.RegularizedLeastSquares
using Test
using LinearAlgebra
using Random
using ImageUtils

# Metal has no Float64, so runtests.jl cannot run here. This runs the ComplexF32 sampled multi-coil
# Cartesian CGNR reconstruction with arrayType = MtlArray and a Float32 regularization weight.

# Load the reconstruction tests without the loop over arrayTypes at the end of the file
include(ex -> Meta.isexpr(ex, :for) ? nothing : ex, joinpath(@__DIR__(), "..", "testReconstruction.jl"))

@testset "MRIReco on Metal" begin
  @testset "Multi Coil CGNR F32" testCSRecoMultiCoilCGNR(type = ComplexF32, arrayType = MtlArray, lambda = 1f-3)
end
