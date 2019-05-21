export CartesianTrajectory3D, cartesian3dNodes, cartesian3dDensity

"""
    CartesianTrajectory3D(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , numSlices=1
                  , kargs...)

returns a 3d cartesian trajectory.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`numSlices=1`)                 - number of slices
"""
function CartesianTrajectory3D(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , numSlices=1
                  , kargs...)
  nodes = cartesian3dNodes(numProfiles, numSamplingPerProfile; numSlices=numSlices)
  times = readoutTimes(numProfiles, numSamplingPerProfile, numSlices; TE=TE, AQ=AQ)
  return  Trajectory("Cartesian3D", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, numSlices, true, false)
end

function cartesian3dNodes(numProfiles, numSamplingPerProfile
                          ; numSlices=1, kargs...)
  nodes = zeros(3,numSamplingPerProfile, numProfiles, numSlices)
  posX = collect( -ceil(Int64, (numSamplingPerProfile-1)/2.):floor(Int64, (numSamplingPerProfile-1)/2.) ) / numSamplingPerProfile
  posY = collect( -ceil(Int64, (numProfiles-1)/2.):floor(Int64, (numProfiles-1)/2.) ) / numProfiles
  posZ = collect( -ceil(Int64, (numSlices-1)/2.):floor(Int64, (numSlices-1)/2.) ) / numSlices

  for j = 1:numSlices
    for l = 1:numProfiles
      for k = 1:numSamplingPerProfile
        nodes[1,k,l,j] = posX[k]
        nodes[2,k,l,j] = posY[l]
        nodes[3,k,l,j] = posZ[j]
      end
    end
  end
  return reshape(nodes, 3, numSamplingPerProfile*numProfiles*numSlices)
end

function cartesian3dDensity(numSamplingPerProfile::Int64, numProfiles::Int64, numSlices::Int64)
  density = zeros(numSamplingPerProfile, numProfiles, numSlices)
  for j = 1:numSlices
    for l = 1:numProfiles
      for k = 1:numSamplingPerProfile
        density[k,l,j] = 1
      end
    end
  end
  return vec(density)
end
