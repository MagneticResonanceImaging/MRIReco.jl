module MRIBase

using Graphics: @mustimplement
using AbstractNFFTs
using NFFTTools # for density compensation weights in trajectory
using SparsityOperators

include("Trajectories/Trajectories.jl")
include("Sequences/Sequence.jl")
include("Datatypes/Datatypes.jl")

end # module
