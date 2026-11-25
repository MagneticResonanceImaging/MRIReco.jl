using Metal
using Test, MRIBase, MRIOperators, MRISimulation, MRIOperators.NFFT, MRIOperators.NFFT.FFTW
using NonuniformFFTs
using LinearAlgebra, MRIOperators.LinearOperatorCollection
using MRIOperators.OhMyThreads
using JLArrays

# Metal has no Float64, so runtests.jl cannot run here. This runs the ComplexF32 Adapt tests:
# SensitivityOp, SubspaceOp and their products moved between Array and MtlArray, compared with the source.
arrayTypes = [MtlArray]
adaptEltypes = (ComplexF32,)

@testset "MRIOperators on Metal" begin
  include(joinpath(@__DIR__(), "..", "testAdapt.jl"))
end
