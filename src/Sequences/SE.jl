export SESequence, echoAmplitudes, flipAngles, numEchoes, TE

mutable struct SESequence <:AbstractSequence
  traj :: AbstractTrajectory    # readout trajectory
end

#
# echo amplitude due to T2-relaxation
#
echoAmplitudes(seq::SESequence, R1::Float64, R2::Float64) = [exp(-seq.traj.TE*R2)]

flipAngles(seq::SESequence) = [90.]

numEchoes(seq::SESequence) = 1

TE(seq::SESequence) = seq.traj.TE

string(seq::SESequence) = "SE"
