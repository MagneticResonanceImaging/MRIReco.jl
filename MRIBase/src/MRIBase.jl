module MRIBase

using Graphics: @mustimplement
using NFFT, NFFTTools # for density compensation weights in trajectory

include("Trajectories/Trajectories.jl")
include("Sequences/Sequence.jl")
include("Datatypes/Datatypes.jl")

end # module
