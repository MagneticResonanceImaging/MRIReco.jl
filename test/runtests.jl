import ImageMagick # not good ...
using MRIReco
using Test
using LinearAlgebra

include("testTrajectories.jl")
include("testIO.jl")
include("testSimulation.jl")
include("testOperators.jl")
include("testReconstruction.jl")

testSimulation()
testOperators()
testReco()

include("testISMRMRD.jl")
