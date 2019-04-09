export RandomCartesianTrajectory, randomCartesianNodes, randomCartesianDensity

function RandomCartesianTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , redFac=1.0
                  , patFunc="poisson"
                  , kargs...)
  nodes = randomCartesianNodes(numProfiles, numSamplingPerProfile; redFac=redFac,patFunc=patFunc)
  times = readoutTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("RandomCartesian", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, true, false)
end

function randomCartesianNodes(numProfiles, numSamplingPerProfile;
                          redFac=1.0, patFunc="poisson", kargs...)
  nodes = zeros(2,numSamplingPerProfile, numProfiles)
  posX = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  posY = collect((0:numProfiles-1)/numProfiles .- 0.5)

  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      nodes[1,k,l] = posX[k]
      nodes[2,k,l] = posY[l]
    end
  end

  nodes = reshape(nodes, 2, numSamplingPerProfile*numProfiles)
  # idx = shuffle_vector(collect(1:numSamplingPerProfile*numProfiles))
  idx = shuffle_vector(collect(1:numSamplingPerProfile*numProfiles), redFac, patFunc)

  return nodes[:,idx]
end

function randomCartesianDensity(numProfiles::Int64, numSamplingPerProfile::Int64)
  density = zeros(numSamplingPerProfile, numProfiles)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      density[k,l] = 1
    end
  end
  return vec(density)
end
