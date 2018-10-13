export StackOfStarsTrajectory

mutable struct StackOfStarsTrajectory <: StackedTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  numSlices::Int
  TE::Float64 # echo time in ms
  AQ::Float64 # time for each spiral arm in ms
  angleOffset::Symbol
end

StackOfStarsTrajectory(numProfiles, numSamplingPerProfile; numSlices=1, TE=0.0, AQ=1.0, angleOffset = :equispaced,kargs...) =
   StackOfStarsTrajectory(numProfiles,numSamplingPerProfile,numSlices,TE,AQ,angleOffset)

string(tr::StackOfStarsTrajectory) = "StackOfStars"


function kspaceNodes(tr::StackOfStarsTrajectory)

  angles = collect(pi.*(0:tr.numProfiles-1)/tr.numProfiles)
  nodes = zeros(3,tr.numSamplingPerProfile, tr.numProfiles, tr.numSlices)
  if tr.angleOffset == :golden
    angles = [i*getGoldenAngleRad()  for i=0:tr.numProfiles-1 ]
  elseif tr.angleOffset == :random
    angles = sort(pi .* rand(tr.numProfiles))
  elseif tr.angleOffset == :equispaced
    angles = collect(pi.*(0:tr.numProfiles-1)/tr.numProfiles)
  end

  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  posZ = collect((0:tr.numSlices-1)/tr.numSlices .- 0.5)
  for j = 1:tr.numSlices
    for l = 1:tr.numProfiles
      for k = 1:tr.numSamplingPerProfile
        nodes[1,k,l,j] = (-1)^l * pos[k]*cos(angles[l])
        nodes[2,k,l,j] = (-1)^l * pos[k]*sin(angles[l])
        nodes[3,k,l,j] = posZ[j]
      end
    end
  end
  return reshape(nodes, 3, tr.numSamplingPerProfile*tr.numProfiles*tr.numSlices)
end

function kspaceDensity(tr::StackOfStarsTrajectory)
  density = zeros(tr.numSamplingPerProfile, tr.numProfiles, tr.numSlices)
  pos = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  for j= 1:tr.numSlices
    for l = 1:tr.numProfiles
      for k = 1:tr.numSamplingPerProfile
        density[k,l,j] = abs(pos[k])
      end
    end
  end
  return vec(density)
end
