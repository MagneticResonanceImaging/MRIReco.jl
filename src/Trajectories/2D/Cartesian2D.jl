export SimpleCartesianTrajectory

# FIXME: Why is this called "Simple?"

mutable struct SimpleCartesianTrajectory <: CartesianTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in s
  AQ::Float64 # time for each spiral arm in s
end

# SimpleCartesianTrajectory(numProfiles, numSamplingPerProfile) =
#    SimpleCartesianTrajectory(numProfiles, numSamplingPerProfile, 0.0, 1.0)

SimpleCartesianTrajectory(numProfiles, numSamplingPerProfile; TE=0.0, AQ=1.0e-3, kargs...) =
   SimpleCartesianTrajectory(numProfiles, numSamplingPerProfile, TE, AQ)

string(tr::SimpleCartesianTrajectory) = "Cartesian"

function kspaceNodes(tr::SimpleCartesianTrajectory)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)
  # posX = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  # posY = collect((0:tr.numProfiles-1)/tr.numProfiles .- 0.5)
  posX = collect( -ceil(Int64, (tr.numSamplingPerProfile-1)/2.):floor(Int64, (tr.numSamplingPerProfile-1)/2.) ) / tr.numSamplingPerProfile
  posY = collect( -ceil(Int64, (tr.numProfiles-1)/2.):floor(Int64, (tr.numProfiles-1)/2.) ) / tr.numProfiles

  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      nodes[1,k,l] = posX[k]
      nodes[2,k,l] = posY[l]
    end
  end
  return reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
end

function kspaceDensity(tr::SimpleCartesianTrajectory)
  density = zeros(tr.numSamplingPerProfile, tr.numProfiles)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      density[k,l] = 1
    end
  end
  return vec(density)
end
