using MRIReco
using Test
using LinearAlgebra

include("testSimulation.jl")
include("testOperators.jl")
include("testReconstruction.jl")

testSimulation()
testOperators()
testReco()
