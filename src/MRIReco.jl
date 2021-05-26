module MRIReco

using Reexport
using FileIO
using ProgressMeter
@reexport using ImageUtils
using NFFT
using Distributions
using LinearOperators
@reexport using RegularizedLeastSquares
using StatsBase
using LinearAlgebra
using Random
@reexport using FFTW
using Distributed
using Graphics: @mustimplement
using Wavelets
using HDF5
using LightXML
using NIfTI
using ThreadPools
@reexport using Unitful
import LowRankApprox.psvd

const Trafo = Union{AbstractMatrix, AbstractLinearOperator, Nothing}

include("Trajectories/Trajectories.jl")
include("Sequences/Sequence.jl")
include("Datatypes/Datatypes.jl")
include("Tools/Tools.jl")
include("Operators/Operators.jl")
include("Simulation/Simulation.jl")
include("Reconstruction/Reconstruction.jl")
include("IO/IO.jl")
include("Sampling/Sampling.jl")

function __init__()
  if Threads.nthreads() > 1
    BLAS.set_num_threads(1)
    FFTW.set_num_threads(1)
  elseif Sys.iswindows()
    BLAS.set_num_threads(1) # see https://github.com/JuliaLang/julia/issues/36976
  end
end

end # module
