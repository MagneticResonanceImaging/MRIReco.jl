using Test, MRIBase, MRIOperators, MRISimulation, MRIOperators.NFFT, MRIOperators.NFFT.FFTW
using LinearAlgebra, MRIOperators.LinearOperatorCollection
using JLArrays


areTypesDefined = @isdefined arrayTypes
arrayTypes = areTypesDefined ? arrayTypes : [Array, JLArray]

@testset "MRIOperators" begin
  include("testOperators.jl")
end