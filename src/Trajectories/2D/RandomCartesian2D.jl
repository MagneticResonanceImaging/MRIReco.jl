export RandomCartesianTrajectory

mutable struct RandomCartesianTrajectory <: CartesianTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in s
  AQ::Float64 # time for each spiral arm in s
  redFac::Float64
  patFunc::AbstractString
end

RandomCartesianTrajectory(numProfiles, numSamplingPerProfile; TE=0.0, AQ=1.0, redFac=1.0, patFunc="poisson", kargs...) =
  RandomCartesianTrajectory(numProfiles, numSamplingPerProfile, TE, AQ, redFac, patFunc)

string(tr::RandomCartesianTrajectory) = "RandomCartesian"

function kspaceNodes(tr::RandomCartesianTrajectory)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)
  posX = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  posY = collect((0:tr.numProfiles-1)/tr.numProfiles .- 0.5)

  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      nodes[1,k,l] = posX[k]
      nodes[2,k,l] = posY[l]
    end
  end

  nodes = reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
  # idx = shuffle_vector(collect(1:tr.numSamplingPerProfile*tr.numProfiles))
  idx = shuffle_vector(collect(1:tr.numSamplingPerProfile*tr.numProfiles), tr.redFac, tr.patFunc)

  return nodes[:,idx]
end

function kspaceDensity(tr::RandomCartesianTrajectory)
  density = zeros(tr.numSamplingPerProfile, tr.numProfiles)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      density[k,l] = 1
    end
  end
  return vec(density)
end
