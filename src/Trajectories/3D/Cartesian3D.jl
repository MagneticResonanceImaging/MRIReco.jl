export CartesianTrajectory3D

mutable struct CartesianTrajectory3D <: StackedTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  numSlices::Int
  TE::Float64 # echo time in s
  AQ::Float64 # time for each spiral arm in s
end

CartesianTrajectory3D(numProfiles, numSamplingPerProfile; numSlices=1, TE=0.0, AQ=1.0e-3, kargs...) =
   CartesianTrajectory3D(numProfiles, numSamplingPerProfile, numSlices, TE, AQ)

string(tr::CartesianTrajectory3D) = "Cartesian3D"


function kspaceNodes(tr::CartesianTrajectory3D)
  nodes = zeros(3,tr.numSamplingPerProfile, tr.numProfiles, tr.numSlices)
  # posX = collect((0:tr.numSamplingPerProfile-1)/tr.numSamplingPerProfile .- 0.5)
  # posY = collect((0:tr.numProfiles-1)/tr.numProfiles .- 0.5)
  # posZ = collect((0:tr.numSlices-1)/tr.numSlices .- 0.5)
  posX = collect( -ceil(Int64, (tr.numSamplingPerProfile-1)/2.):floor(Int64, (tr.numSamplingPerProfile-1)/2.) ) / tr.numSamplingPerProfile
  posY = collect( -ceil(Int64, (tr.numProfiles-1)/2.):floor(Int64, (tr.numProfiles-1)/2.) ) / tr.numProfiles
  posZ = collect( -ceil(Int64, (tr.numSlices-1)/2.):floor(Int64, (tr.numSlices-1)/2.) ) / tr.numSlices

  for j = 1:tr.numSlices
    for l = 1:tr.numProfiles
      for k = 1:tr.numSamplingPerProfile
        nodes[1,k,l,j] = posX[k]
        nodes[2,k,l,j] = posY[l]
        nodes[3,k,l,j] = posZ[j]
      end
    end
  end
  return reshape(nodes, 3, tr.numSamplingPerProfile*tr.numProfiles*tr.numSlices)
end

function kspaceDensity(tr::CartesianTrajectory3D)
  density = zeros(tr.numSamplingPerProfile, tr.numProfiles, tr.numSlices)
  for j = 1:tr.numSlices
    for l = 1:tr.numProfiles
      for k = 1:tr.numSamplingPerProfile
        density[k,l,j] = 1
      end
    end
  end
  return vec(density)
end
