export CustomTrajectory

mutable struct CustomTrajectory <: Abstract2DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in s
  AQ::Float64 # time for each spiral arm in s
  nodes::Vector{Float64}
  readoutTimes::Vector{Float64}
  isCirc::Bool
end

function CustomTrajectory(numProfiles, numSamplingPerProfile, nodes; readoutTimes=nothing, 
				       TE=0.0, AQ=1.0e-3, isCircular=false, kargs...)
  if readoutTimes != nothing
    times = readoutTimes
  else
    times = zeros(numSamplingPerProfile, numProfiles)
    for l = 1:numProfiles
      for k = 1:numSamplingPerProfile
        times[k,l] = TE + AQ*(k-1)/numSamplingPerProfile
      end
    end
  end
  CustomTrajectory(numProfiles, numSamplingPerProfile, TE, AQ, nodes, vec(times), isCircular)
end

string(tr::CustomTrajectory) = "Custom"

kspaceNodes(tr::CustomTrajectory) = reshape(tr.nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)

readoutTimes(tr::CustomTrajectory) = tr.readoutTimes

isCircular(tr::CustomTrajectory) = tr.isCirc
