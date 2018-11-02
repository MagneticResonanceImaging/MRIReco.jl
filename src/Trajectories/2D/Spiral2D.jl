export SpiralTrajectory

# TODO: In this file voronoi based density calculation is done. This is,
# however not restricted to Spirals and should be moved out.

mutable struct SpiralTrajectory <: Abstract2DTrajectory
  numProfiles::Int64
  numSamplingPerProfile::Int64
  TE::Float64 # echo time in ms
  AQ::Float64 # time for each spiral arm in ms
  windings::Real
  angleOffset::String
end

function SpiralTrajectory(numProfiles::Int64
                            , numSamplingPerProfile::Int64
                            ; TE=0.0
                            , AQ=1.0e-3
                            , windings= 6.25
                            , angleOffset= "equispaced"
                            , randomSampleOnProfile=false
                            , kargs...)

    SpiralTrajectory(numProfiles
                    ,numSamplingPerProfile
                    , TE
                    , AQ
                    , windings
                    , angleOffset)
end

string(tr::SpiralTrajectory) = "Spiral"

function kspaceNodes(tr::SpiralTrajectory)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)
  if tr.angleOffset == "golden"
      angles = [i*(3-sqrt(5))/2  for i=0:tr.numProfiles-1 ]
  elseif tr.angleOffset == "random"
      angles = sort(rand(tr.numProfiles))
  elseif tr.angleOffset == "equispaced"
      angles = collect((0:tr.numProfiles-1)/tr.numProfiles)
  end

  A = 0.5 # Maximum radius is 0.5
  w = tr.windings # 8/64 * 50 = 6.25 which means we have 6.25 turns aka windings
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      # For the special case that τ(t)=t, the amount of time spent for each
      # winding is constant, regardless of whether the acquired
      # winding is near the center or in the outer part of the spiral. I
      # other words, the readout gradients reach their maximum
      # performance at the end of the acquisition.
      t = sqrt((k-1)/(tr.numSamplingPerProfile-1)) #
      nodes[1,k,l] = A*t*cos(2*pi*( w*t + angles[l] ))
      nodes[2,k,l] = A*t*sin(2*pi*( w*t + angles[l] ))
    end
  end
  return reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
end

"""
Calculating the density using VoronoiCells
"""
function kspaceDensity(tr::SpiralTrajectory)
  error("Not implemented!")
end

isCircular(tr::SpiralTrajectory) = true
