export RandomCartesianTrajectory, randomCartesianNodes, randomCartesianDensity

"""
    CartesianTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , kmin=(-0.5,-0.5)
                  , kmax=(0.5,0.5)
                  , kargs...)

returns a 2d cartesian trajectory with random subsampling

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`redFac=1.0)`)                 - subsampling Factor
* (`patFunc="poisson"`)           - name of sampling pattern to use
"""
function RandomCartesianTrajectory(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , redFac=1.0
                  , patFunc="poisson"
                  , kargs...) where T
  nodes = randomCartesianNodes(T, numProfiles, numSamplingPerProfile; redFac=redFac,patFunc=patFunc)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("RandomCartesian", nodes, times, T(TE), T(AQ), numProfiles, numSamplingPerProfile, 1, true, false)
end

function randomCartesianNodes(::Type{T}, numProfiles, numSamplingPerProfile;
                          redFac=1.0, patFunc="poisson", kargs...) where T
  nodes = zeros(T, 2,numSamplingPerProfile, numProfiles)
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

function randomCartesianDensity(::Type{T}, numProfiles::Int64, numSamplingPerProfile::Int64) where T
  density = zeros(T, numSamplingPerProfile, numProfiles)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      density[k,l] = 1
    end
  end
  return vec(density)
end
