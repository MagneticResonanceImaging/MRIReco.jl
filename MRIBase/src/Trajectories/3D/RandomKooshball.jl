export RandomKooshballTrajectory, randomKooshballNodes, randomKooshballDensity

"""
    KooshballTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , kargs...)

returns a 3d kooshball trajectory with randomized order of the profiles

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
"""
function RandomKooshballTrajectory(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , kargs...) where T
  nodes = randomKooshballNodes(T, numProfiles, numSamplingPerProfile; numSlices=1)
  times = kooshballTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("RandomKooshball", nodes, times, T(TE), T(AQ), numProfiles, numSamplingPerProfile, 1, false, true)
end

function randomKooshballNodes(::Type{T}, numProfiles, numSamplingPerProfile; kargs...) where T
  nodes = zeros(T, 3, numSamplingPerProfile, numProfiles)
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  profileList = shuffle_vector([i for i=1:numProfiles])

  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      nodes[1,k,l] = (-1)^l * pos[k]*xax(profileList[l],2*numProfiles)
      nodes[2,k,l] = (-1)^l * pos[k]*yax(profileList[l],2*numProfiles)
      nodes[3,k,l] = (-1)^l * pos[k]*zax(profileList[l],2*numProfiles)
    end
  end
  return reshape(nodes, 3, numSamplingPerProfile*numProfiles)
end

function randomKooshballDensity(::Type{T}, numProfiles::Int64, numSamplingPerProfile::Int64) where T
  density = zeros(T, numSamplingPerProfile, numProfiles)
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      density[k,l] = abs(pos[k])
    end
  end
  return vec(density)
end
