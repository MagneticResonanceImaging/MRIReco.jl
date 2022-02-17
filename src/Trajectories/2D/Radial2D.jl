export RadialTrajectory, radialNodes, radialDensity

"""
    RadialTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , angleOffset= :equispaced
                  , kargs...)

returns a 2d radial trajectory.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`angleOffset= :equispaced`)    - spacing of profile angles (`:equispaced` sampling, `:golden` angle sampling or `:random` sampling)
"""
function RadialTrajectory(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE::T=0.0
                  , AQ::T=1.e-3
                  , angleOffset= :equispaced
                  , kargs...) where T
  nodes = radialNodes(T, numProfiles, numSamplingPerProfile; angleOffset=angleOffset)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("Radial", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

function radialNodes(::Type{T}, numProfiles::Int64
                          , numSamplingPerProfile::Int64
                          ; angleOffset= :equispaced
                          , kargs...) where T
  nodes = zeros(T, 2, numSamplingPerProfile, numProfiles)
  if angleOffset == :golden
    angles = [i*getGoldenAngleRad()  for i=0:numProfiles-1 ]
  elseif angleOffset == :random
    angles = sort(pi .* rand(numProfiles))
  elseif angleOffset == :equispaced
    angles = collect(pi.*(0:numProfiles-1)/numProfiles)
  end


  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      nodes[1,k,l] = (-1)^l * pos[k]*cos(angles[l])
      nodes[2,k,l] = (-1)^l * pos[k]*sin(angles[l])
    end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end

function radialDensity(::Type{T}, numProfiles::Int64, numSamplingPerProfile::Int64) where T
  density = zeros(T, numSamplingPerProfile, numProfiles)
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      density[k,l] = abs(pos[k])
    end
  end
  return vec(density)
end
