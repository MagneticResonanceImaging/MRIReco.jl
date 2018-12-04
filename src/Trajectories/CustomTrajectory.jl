export CustomTrajectory, CustomTrajectory2D, CustomTrajectory3D

mutable struct CustomTrajectory2D <: Abstract2DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in s
  AQ::Float64 # time for each spiral arm in s
  nodes::Vector{Float64}
  readoutTimes::Vector{Float64}
  isCirc::Int64
end

mutable struct CustomTrajectory3D <: Abstract3DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  numSlices::Int
  TE::Float64 # echo time in s
  AQ::Float64 # time for each spiral arm in s
  nodes::Vector{Float64}
  readoutTimes::Vector{Float64}
  isCirc::Int64
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
  CustomTrajectory2D(numProfiles, numSamplingPerProfile, TE, AQ, nodes, vec(times), Int64(isCircular))
end

function CustomTrajectory(numProfiles, numSamplingPerProfile, numSlices, nodes; readoutTimes=nothing,
				       TE=0.0, AQ=1.0e-3, isCircular=false, kargs...)
  if readoutTimes != nothing
    times = readoutTimes
  else
    times = zeros(numSamplingPerProfile, numProfiles, numSlices)
    for m = 1:numSlices
      for l = 1:numProfiles
        for k = 1:numSamplingPerProfile
          times[k,l,m] = TE + AQ*(k-1)/numSamplingPerProfile
        end
      end
    end
  end
  CustomTrajectory3D(numProfiles, numSamplingPerProfile, numSlices, TE, AQ, nodes, vec(times), Int64(isCircular))
end

string(tr::CustomTrajectory2D) = "Custom2D"
string(tr::CustomTrajectory3D) = "Custom3D"

kspaceNodes(tr::CustomTrajectory2D) = reshape(tr.nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
kspaceNodes(tr::CustomTrajectory3D) = reshape(tr.nodes, 3, tr.numSamplingPerProfile*tr.numProfiles)

readoutTimes(tr::T) where T <: Union{CustomTrajectory2D,CustomTrajectory3D} = tr.readoutTimes

isCircular(tr::T) where T <: Union{CustomTrajectory2D,CustomTrajectory3D} = (tr.isCirc==1)
