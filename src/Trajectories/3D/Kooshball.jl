export KooshballTrajectory

mutable struct KooshballTrajectory <: Real3DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in ms
  AQ::Float64 # time for each spiral arm in ms
end

KooshballTrajectory(numProfiles, numSamplingPerProfile; TE=0.0, AQ=1.0, kargs...) =
    KooshballTrajectory(numProfiles, numSamplingPerProfile, TE, AQ)

string(tr::KooshballTrajectory) = "Kooshball"

function kspaceNodes(tr::KooshballTrajectory)
  nodes = zeros(3,tr.numSamplingPerProfile, tr.numProfiles)
  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      nodes[1,k,l] = (-1)^l * pos[k]*xax(l,2*tr.numProfiles)
      nodes[2,k,l] = (-1)^l * pos[k]*yax(l,2*tr.numProfiles)
      nodes[3,k,l] = (-1)^l * pos[k]*zax(l,2*tr.numProfiles)
    end
  end
  return reshape(nodes, 3, tr.numSamplingPerProfile*tr.numProfiles)
end

function kspaceDensity(tr::KooshballTrajectory)
  density = zeros(tr.numSamplingPerProfile, tr.numProfiles)
  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      density[k,l] = abs(pos[k])
    end
  end
  return vec(density)
end
