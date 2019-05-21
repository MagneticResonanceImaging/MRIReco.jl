export KooshballTrajectory, kooshballNodes, kooshballDensity, kooshBallTimes

"""
    KooshballTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , kargs...)

returns a 3d kooshball trajectory.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
"""
function KooshballTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , kargs...)
  nodes = kooshballNodes(numProfiles, numSamplingPerProfile; numSlices=numSlices)
  times = kooshballTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("Kooshball", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

function kooshballNodes(numProfiles::Int64, numSamplingPerProfile::Int64; kargs...)
  nodes = zeros(3,numSamplingPerProfile, numProfiles)
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      nodes[1,k,l] = (-1)^l * pos[k]*xax(l,2*numProfiles)
      nodes[2,k,l] = (-1)^l * pos[k]*yax(l,2*numProfiles)
      nodes[3,k,l] = (-1)^l * pos[k]*zax(l,2*numProfiles)
    end
  end
  return reshape(nodes, 3, numSamplingPerProfile*numProfiles)
end

function kooshballTimes(numProfiles::Int64, numSamplingPerProfile::Int64; TE=0.0, AQ=1.e-3, kargs...)
  times = zeros(numSamplingPerProfile, numProfiles)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      times[k,l] = (TE + AQ*(k-1)/numSamplingPerProfile )
    end
  end
  return return vec(times)
end

function kooshballDensity(numSamplingPerProfile::Int64, numProfiles::Int64, numSlices::Int64)
  density = zeros(numSamplingPerProfile, numProfiles)
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      density[k,l] = abs(pos[k])
    end
  end
  return vec(density)
end
