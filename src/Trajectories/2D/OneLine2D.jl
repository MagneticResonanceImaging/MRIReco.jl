export OneLine2Trajectory, oneLine2dNodes, oneLine2dDensity

function OneLine2dTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , angle::Float64=0.0
                  , kargs...)
  nodes = oneLine2dNodes(numProfiles, numSamplingPerProfile; angle=angle)
  times = readoutTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("OneLine", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

function oneLine2dNodes(numProfiles, numSamplingPerProfile; angle::Float64=0.0, kargs...)
  nodes = zeros(2,numSamplingPerProfile, numProfiles)
  angles = zeros(Float64,numProfiles)
  angles[:] .= angle
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      nodes[1,k,l] = (-1)^l * pos[k]*cos(angles[l])
      nodes[2,k,l] = (-1)^l * pos[k]*sin(angles[l])
    end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end

function oneLine2dDensity(numSamplingPerProfile::Int64, numProfiles::Int64)
  density = zeros(tr.numSamplingPerProfile, numProfiles)
  pos = collect((0:numSamplingPerProfile-1)/numSamplingPerProfile .- 0.5)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      density[k,l] = abs(pos[k])
    end
  end
  return vec(density)
end
