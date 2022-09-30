export SpiralTrajectoryDualDensity

"""
    SpiralTrajectoryDualDensity(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , windings::Real= 6.25
                  , angleOffset= :equispaced
                  , densityFactor=2.0
                  , kargs...)

returns a 2d spiral trajectory.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`windings::Real= 6.25`)        - number of windings of the spiral profiles
* (`angleOffset= :equispaced`)    - spacing of profile angles (`:equispaced` sampling, `:golden` angle sampling or `:random` sampling)
"""
function SpiralTrajectoryDualDens(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , windings::Real= 6.25
                  , kmax::Real=0.5
                  , angleOffset::String="equispaced"
                  , densityFactor::Float64=2.0
                  , kargs...) where T
  nodes = spiralNodesDualDens(T, numProfiles, numSamplingPerProfile; windings=windings,
              angleOffset=angleOffset, densityFactor=densityFactor, kmax=kmax)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("SpiralDualDens", nodes, times, T(TE), T(AQ), numProfiles, numSamplingPerProfile, 1, false, true)
end

function spiralNodesDualDens(::Type{T}, numProfiles::Int64
                            , numSamplingPerProfile::Int64
                            ; windings::Real= 6.25
                            , kmax::Real=0.5
                            , angleOffset::String="equispaced"
                            , densityFactor::Float64
                            , kargs...) where T

  nodes = zeros(T, 2, numSamplingPerProfile, numProfiles)
  if angleOffset == "golden"
      angles = [i*(3-sqrt(5))/2  for i=0:numProfiles-1 ]
  elseif angleOffset == "random"
      angles = sort(rand(numProfiles))
  elseif angleOffset == "equispaced"
      angles = collect((0:numProfiles-1)/numProfiles)
  end

  w = windings
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile

      t = sqrt((k-1)/(numSamplingPerProfile-1))

      if t < 1/2
        tau = t/densityFactor
      else
        tau = 0.5/densityFactor + (t-0.5)*2*(1.0-0.5/densityFactor)
      end

      nodes[1,k,l] = kmax*tau*cos(2*pi*( w*t + angles[l] ))
      nodes[2,k,l] = kmax*tau*sin(2*pi*( w*t + angles[l] ))
    end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end
