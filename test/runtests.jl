using MRIReco
using Test
using LinearAlgebra

include("testIO.jl")
include("testSimulation.jl")
include("testOperators.jl")
include("testReconstruction.jl")

testSimulation()
testOperators()
testReco()

include("testISMRMRD.jl")
