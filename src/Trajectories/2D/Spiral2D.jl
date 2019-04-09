export SpiralTrajectory, spiralNodes, spiralDensity

function SpiralTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , windings::Real= 6.25
                  , angleOffset::String="equispaced"
                  , kargs...)
  nodes = spiralNodes(numProfiles, numSamplingPerProfile; windings=windings, angleOffset=angleOffset)
  times = readoutTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("Spiral", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

function spiralNodes(numProfiles::Int64
                            , numSamplingPerProfile::Int64
                            ; windings::Real= 6.25
                            , angleOffset::String="equispaced"
                            , kargs...)

  nodes = zeros(2,numSamplingPerProfile, numProfiles)
  if angleOffset == "golden"
      angles = [i*(3-sqrt(5))/2  for i=0:numProfiles-1 ]
  elseif angleOffset == "random"
      angles = sort(rand(numProfiles))
  elseif angleOffset == "equispaced"
      angles = collect((0:numProfiles-1)/numProfiles)
  end

  A = 0.5 # Maximum radius is 0.5
  w = windings # 8/64 * 50 = 6.25 which means we have 6.25 turns aka windings
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      # For the special case that τ(t)=t, the amount of time spent for each
      # winding is constant, regardless of whether the acquired
      # winding is near the center or in the outer part of the spiral. I
      # other words, the readout gradients reach their maximum
      # performance at the end of the acquisition.
      t = sqrt((k-1)/(numSamplingPerProfile-1)) #
      nodes[1,k,l] = A*t*cos(2*pi*( w*t + angles[l] ))
      nodes[2,k,l] = A*t*sin(2*pi*( w*t + angles[l] ))
    end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end

"""
Calculating the density using VoronoiCells
"""
function spiralDensity(numProfiles::Int64, numSamplingPerProfile::Int64)
  @error "Not implemented!"
end
