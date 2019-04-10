export AbstractSequence, AbstractMultiEchoSequence, echoAmplitudes, numEchoes, echoTimes, sequence, sequenceInfo

abstract type AbstractSequence end

# should contain a method to calculate echo amplitudes
# trajectories should be specified for each echo
abstract type AbstractMultiEchoSequence <: AbstractSequence end

# define interface for a sequence
@mustimplement trajectory(seq::AbstractSequence,n::Int64)
@mustimplement echoAmplitudes(seq::AbstractSequence,R1::Float64,R2::Float64)
@mustimplement numEchoes(seq::AbstractSequence)
@mustimplement encoding(seq::AbstractSequence)
@mustimplement string(seq::AbstractSequence)

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

function sequenceInfo(seq::AbstractSequence)
  param = Dict{Symbol,Any}()
  param[:name] = string(seq)
  for field in fieldnames(typeof(seq))
    param[field] = getfield(seq,field)
  end
  return param
end
