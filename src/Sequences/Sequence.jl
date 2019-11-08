export AbstractSequence, AbstractMultiEchoSequence, echoAmplitudes, numEchoes, echoTimes, sequence, sequenceInfo

abstract type AbstractSequence end

# define interface for a sequence
@mustimplement echoAmplitudes(seq::AbstractSequence,R1::Float64,R2::Float64)
@mustimplement numContrasts(seq::AbstractSequence)
@mustimplement string(seq::AbstractSequence)

include("BlochSimulation.jl")
include("MultiEcho.jl")

function sequence(seqName::AbstractString, numProfiles, numSamplingPerProfile; kargs...)
  if seqName == "ME"
    return MESequence(; kargs...)
  else
    error("The sequence $seqName is not yet supported!")
  end

  return MESequence(; kargs...)
end

function sequenceInfo(seq::AbstractSequence)
  param = Dict{Symbol,Any}()
  param[:name] = string(seq)
  for field in fieldnames(typeof(seq))
    param[field] = getfield(seq,field)
  end
  return param
end
