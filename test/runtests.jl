using MRIReco
using Test
using LinearAlgebra
using HTTP
using Random
using StatsBase

using Scratch
using LazyArtifacts

const datadir = joinpath(artifact"MRIRecoTestData")
@info "The test data is located at $datadir."

const tmpdir  = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpdir."


include("testTrajectories.jl")
include("testIO.jl")
include("testSimulation.jl")
include("testOperators.jl")
include("testReconstruction.jl")
include("testBrukerFile.jl")
include("testISMRMRD.jl")
include("testCoilsens.jl")
