using Test, MRIBase, MRIOperators, MRISimulation, NFFT.FFTW
using LinearAlgebra, LinearOperatorCollection
using JLArrays

arrayTypes = [Array] # , JLArray]

@testset "MRIOperators" begin
  include("testOperators.jl")
end