@parameter struct EmptyForTests <: AbstractMRIRecoParameters end
@reconstruction struct PlacerHolderAlgo <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::EmptyForTests
end

include("testMRIRecoContext.jl")
include("Parameters/testParameters.jl")
include("testIterativeContextParameters.jl")
#include("testPlanReconstructions.jl")