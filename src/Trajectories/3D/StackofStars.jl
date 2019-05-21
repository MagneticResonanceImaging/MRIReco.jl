export StackOfStarsTrajectory, stackOfStarsNodes, stackOfStarsDensity

"""
    StackOfStarsTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , numSlices=1
                  , angleOffset= :equispaced
                  , kargs...)

returns a 2d radial trajectory.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`numSlices=1`)                 - number of slices
* (`angleOffset= :equispaced`)    - spacing of profile angles (`:equispaced` sampling, `:golden` angle sampling or `:random` sampling)
"""
function StackOfStarsTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , numSlices=1
                  , angleOffset = :equispaced
                  , kargs...)
  nodes = cartesian3dNodes(numProfiles, numSamplingPerProfile; numSlices=numSlices, angleOffset = :equispaced)
  times = readoutTimes(numProfiles, numSamplingPerProfile, numSlices; TE=TE, AQ=AQ)
  return  Trajectory("StackOfStars", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, numSlices, false, true)
end

function stackOfStarsNodes(numProfiles, numSamplingPerProfile
              ; numSlices=1, angleOffset = :equispaced, kargs...)

  angles = collect(pi.*(0:numProfiles-1)/numProfiles)
  nodes = zeros(3,numSamplingPerProfile, numProfiles, numSlices)
  if angleOffset == :golden
    angles = [i*getGoldenAngleRad()  for i=0:numProfiles-1 ]
  elseif angleOffset == :random
    angles = sort(pi .* rand(numProfiles))
  elseif angleOffset == :equispaced
    angles = collect(pi.*(0:numProfiles-1)/numProfiles)
  end

  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  posZ = collect((0:numSlices-1)/numSlices .- 0.5)
  for j = 1:numSlices
    for l = 1:numProfiles
      for k = 1:numSamplingPerProfile
        nodes[1,k,l,j] = (-1)^l * pos[k]*cos(angles[l])
        nodes[2,k,l,j] = (-1)^l * pos[k]*sin(angles[l])
        nodes[3,k,l,j] = posZ[j]
      end
    end
  end
  return reshape(nodes, 3, numSamplingPerProfile*numProfiles*numSlices)
end

function stackOfStarsDensity(numSamplingPerProfile::Int64, numProfiles::Int64, numSlices::Int64)
  density = zeros(numSamplingPerProfile, numProfiles, numSlices)
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for j= 1:numSlices
    for l = 1:numProfiles
      for k = 1:numSamplingPerProfile
        density[k,l,j] = abs(pos[k])
      end
    end
  end
  return vec(density)
end
