module MRIBase

using Graphics: @mustimplement
using AbstractNFFTs
using NFFTTools # for density compensation weights in trajectory

include("Trajectories/Trajectories.jl")
include("Sequences/Sequence.jl")
include("Datatypes/Datatypes.jl")

"""
    correctOffset(acq::AcquisitionData,offsetCor=[0,0,0])

    Correct in the k-space the offset along read/phase1/phase2

"""
function correctOffset(acq::AcquisitionData,offsetCor=[0,0,0])
    contrs = size(acq.kdata,1)
    sls = size(acq.kdata,2)
    reps = size(acq.kdata,3)

    shift = offsetCor ./ acq.fov
    shift = shift .* float.(acq.encodingSize)
    # nodes -> -0.5 to 0.5
    for k = 1:reps
        phase_nodes = exp.(-2π * im * acq.traj[k].nodes' * shift)
        [acq.kdata[i,j,k] = acq.kdata[i,j,k] .* phase_nodes for i=1:contrs, j=1:sls]
    end
    return acq
end

end # module
