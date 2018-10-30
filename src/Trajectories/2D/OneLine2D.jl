export OneLine2DTrajectory

mutable struct OneLine2DTrajectory <: Abstract2DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  angle::Float64
  TE::Float64 # echo time in ms
  AQ::Float64 # time for each spiral arm in ms
end

OneLine2DTrajectory(numProfiles, numSamplingPerProfile;angle=0.0, TE=0.0, AQ=1.0, kargs...) =
   OneLine2DTrajectory(numProfiles, numSamplingPerProfile, angle,TE, AQ)

string(tr::OneLine2DTrajectory) = "OneLine"

function kspaceNodes(tr::OneLine2DTrajectory)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)
  angles = zeros(Float64,tr.numProfiles)
  angles[:] .= tr.angle
  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      nodes[1,k,l] = (-1)^l * pos[k]*cos(angles[l])
      nodes[2,k,l] = (-1)^l * pos[k]*sin(angles[l])
    end
  end
  return reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
end

function kspaceDensity(tr::OneLine2DTrajectory)
  density = zeros(tr.numSamplingPerProfile, tr.numProfiles)
  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      density[k,l] = abs(pos[k])
    end
  end
  return vec(density)
end
