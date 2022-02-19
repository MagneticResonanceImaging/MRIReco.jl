export SpiralTrajectory, spiralNodes, spiralDensity

"""
    SpiralTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , windings::Real= 6.25
                  , angleOffset= :equispaced
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
function SpiralTrajectory(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , windings::Real= 6.25
                  , kmax::Real=0.5
                  , angleOffset::String="equispaced"
                  , kargs...) where T
  nodes = spiralNodes(T, numProfiles, numSamplingPerProfile; windings=windings,
                      angleOffset=angleOffset, kmax=kmax)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("Spiral", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

function spiralNodes(::Type{T}, numProfiles::Int64
                            , numSamplingPerProfile::Int64
                            ; windings::Real= 6.25
                            , kmax::Real=0.5
                            , angleOffset::String="equispaced"
                            , kargs...) where T

  nodes = zeros(T, 2, numSamplingPerProfile, numProfiles)
  if angleOffset == "golden"
      angles = [i*(3-sqrt(5))/2  for i=0:numProfiles-1 ]
  elseif angleOffset == "random"
      angles = sort(rand(numProfiles))
  elseif angleOffset == "equispaced"
      angles = collect((0:numProfiles-1)/numProfiles)
  end

  w = windings # 8/64 * 50 = 6.25 which means we have 6.25 turns aka windings
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      # For the special case that τ(t)=t, the amount of time spent for each
      # winding is constant, regardless of whether the acquired
      # winding is near the center or in the outer part of the spiral. I
      # other words, the readout gradients reach their maximum
      # performance at the end of the acquisition.
      t = sqrt((k-1)/(numSamplingPerProfile-1)) #
      nodes[1,k,l] = kmax*t*cos(2*pi*( w*t + angles[l] ))
      nodes[2,k,l] = kmax*t*sin(2*pi*( w*t + angles[l] ))
    end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end

"""
Calculating the density using VoronoiCells
"""
function spiralDensity(::Type{T}, numProfiles::Int64, numSamplingPerProfile::Int64) where T
  @error "Not implemented!"
end
