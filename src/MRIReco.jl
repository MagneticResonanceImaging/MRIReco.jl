module MRIReco

using Reexport

@reexport using MRIBase
@reexport using MRIOperators

using ProgressMeter
using AxisArrays
using RegularizedLeastSquares
using LinearAlgebra
using Random
using FLoops
using Unitful


include("Tools/Tools.jl")
include("Reconstruction/Reconstruction.jl")

function __init__()
  if Threads.nthreads() > 1
    BLAS.set_num_threads(1)
    FFTW.set_num_threads(1)
  elseif Sys.iswindows()
    BLAS.set_num_threads(1) # see https://github.com/JuliaLang/julia/issues/36976
  end
end

#include("precompile.jl")

end # module

