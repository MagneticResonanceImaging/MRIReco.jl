using Test, MRIBase, MRIOperators, MRISimulation, MRIOperators.NFFT, MRIOperators.NFFT.FFTW
using NonuniformFFTs
using LinearAlgebra, MRIOperators.LinearOperatorCollection
using MRIOperators.OhMyThreads
#using JLArrays


areTypesDefined = @isdefined arrayTypes
arrayTypes = areTypesDefined ? arrayTypes : [Array] # , JLArray]

@testset "MRIOperators" begin
  for backend in [NFFT.backend(), NonuniformFFTs.backend()]
    with(nfft_backend => backend) do
      @testset "Operators with $(string(typeof(backend)))" begin
        include("testOperators.jl")
      end
    end
  end
end