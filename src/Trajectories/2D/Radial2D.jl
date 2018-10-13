export RadialTrajectory

mutable struct RadialTrajectory <: Abstract2DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in ms
  AQ::Float64 # time for each spiral arm in ms
  angleOffset::Symbol
end

function RadialTrajectory(numProfiles::Int64
                          , numSamplingPerProfile::Int64
                          ; TE=0.0
                          , AQ=1.0
                          , angleOffset= :equispaced
                          , kargs...)

   RadialTrajectory(numProfiles, numSamplingPerProfile, TE, AQ, angleOffset)
end

string(tr::RadialTrajectory) = "Radial"

function kspaceNodes(tr::RadialTrajectory)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)
  if tr.angleOffset == :golden
    angles = [i*getGoldenAngleRad()  for i=0:tr.numProfiles-1 ]
  elseif tr.angleOffset == :random
    angles = sort(pi .* rand(tr.numProfiles))
  elseif tr.angleOffset == :equispaced
    angles = collect(pi.*(0:tr.numProfiles-1)/tr.numProfiles)
  end


  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      nodes[1,k,l] = (-1)^l * pos[k]*cos(angles[l])
      nodes[2,k,l] = (-1)^l * pos[k]*sin(angles[l])
    end
  end
  return reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
end

function kspaceDensity(tr::RadialTrajectory)
  density = zeros(tr.numSamplingPerProfile, tr.numProfiles)
  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      density[k,l] = abs(pos[k])
    end
  end
  return vec(density)
end
