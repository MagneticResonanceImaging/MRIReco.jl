using Test, LinearAlgebra, MRIFiles, MRIReco, MRICoilSensitivities
using MRIReco.RegularizedLeastSquares
using Scratch
using LazyArtifacts
using ImageUtils

const datadir = joinpath(artifact"MRIRecoTestData")
@info "The test data is located at $datadir."

const tmpdir  = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpdir."

include("testBrukerFile.jl")
include("testISMRMRD.jl")