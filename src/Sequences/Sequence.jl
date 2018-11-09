export AbstractSequence, AbstractMultiEchoSequence, echoAmplitudes, numEchoes, echoTimes, sequence, setTrajectory!

abstract type AbstractSequence end

# should contain a method to calculate echo amplitudes
# trajectories should be specified for each echo
abstract type AbstractMultiEchoSequence <: AbstractSequence end

include("BlochSimulation.jl")
include("FSE.jl")
include("SE.jl")
include("FISP.jl")

function sequence(seqName::AbstractString, trajName::AbstractString, numProfiles, numSamplingPerProfile; kargs...)
  tr = trajectory(trajName, numProfiles, numSamplingPerProfile; kargs...)
  if seqName == "SE"
    return SESequence(tr)
  elseif seqName == "FSE"
    return FSESequence(tr; kargs...)
  elseif seqName == "FISP"
    return FISPSequence(tr; kargs...)
  else
    error("The sequence $seqName is not yet supported!")
  end

  return SESequence(tr; kargs...)
end


# echo Amplitudes for a standard spin echo sequence
# echoAmplitudes(seq::AbstractTrajectory, T1::Float64, T2::Float64) = [exp(-seq.traj.TE/T2)]

trajectory(seq::AbstractSequence, n=1) = seq.traj

trajectory(seq::AbstractMultiEchoSequence, n=1) = seq.traj[n]

function setTrajectory!(seq::AbstractSequence,tr::AbstractTrajectory, n=1)
  seq.traj = tr
end

function setTrajectory!(seq::AbstractMultiEchoSequence,tr::AbstractTrajectory, n=1)
  seq.traj[n] = tr
end

numEchoes(seq::AbstractSequence) = 1

numEchoes(seq::AbstractMultiEchoSequence) = seq.numEchoes

echoTimes(seq::AbstractSequence) = seq.traj.TE

echoTimes(seq::AbstractMultiEchoSequence) = [ i*seq.traj.TE for i=1:numEchoes(seq) ]
