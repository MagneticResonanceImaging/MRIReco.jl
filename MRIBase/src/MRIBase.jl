module MRIBase

using AbstractNFFTs
using NFFTTools # for density compensation weights in trajectory


include("Trajectories/Trajectories.jl")
include("Datatypes/Datatypes.jl")
include("Datatypes/Flags.jl")
end # module
