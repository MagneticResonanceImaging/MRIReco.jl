module MRIBase

using Graphics: @mustimplement
using NFFTTools # for density compensation weights in trajectory

include("Trajectories/Trajectories.jl")
include("Sequences/Sequence.jl")
include("Datatypes/Datatypes.jl")

end # module
