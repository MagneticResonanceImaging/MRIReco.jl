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
function CartesianTrajectory3D(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , numSlices=1
                  , kargs...) where T
  nodes = cartesian3dNodes(T, numProfiles, numSamplingPerProfile; numSlices=numSlices)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile, numSlices; TE=TE, AQ=AQ)
  return  Trajectory("Cartesian3D", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, numSlices, true, false)
end

function cartesian3dNodes(::Type{T}, numProfiles, numSamplingPerProfile
                          ; numSlices=1, kargs...) where T
  nodes = zeros(T, 3, numSamplingPerProfile, numProfiles, numSlices)
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

function cartesianSubsamplingIdx(shape::NTuple{3,Int}, tr::Trajectory)
  if !isCartesian(tr)
    @error "SampledFFTOp can only be applied to Cartesian data"
  end
  
  nx,ny,nz = shape
  numProf, numSamp, numSl = (numProfiles(tr), numSamplingPerProfile(tr), numSlices(tr))
  idxX = range(1,step=floor(Int, nx/numSamp), length=numSamp)
  idxY = range(1,step=floor(Int, ny/numProf), length=numProf)
  idxZ = range(1,step=floor(Int, nz/numSl), length=numSl)

  idx = Vector{Int}(undef,numProf*numSamp*numSl)
  linIdx1 = LinearIndices((nx,ny,nz))
  linIdx2 = LinearIndices((numSamp,numProf,numSl))
  for k=1:numSl, j=1:numProf, i=1:numSamp
    idx[linIdx2[i,j,k]] = linIdx1[idxX[i],idxY[j],idxZ[k]]
  end

  return idx
end

function isUndersampledCartTrajectory(shape::NTuple{3,Int}, tr::Trajectory)
  return shape[1]!=numSamplingPerProfile(tr) || shape[2]!=numProfiles(tr) || shape[3]!=numSlices(tr)
end

function cartesian3dDensity(::Type{T}, numSamplingPerProfile::Int64, numProfiles::Int64, numSlices::Int64) where T
  density = zeros(T, numSamplingPerProfile, numProfiles, numSlices)
  for j = 1:numSlices
    for l = 1:numProfiles
      for k = 1:numSamplingPerProfile
        density[k,l,j] = 1
      end
    end
  end
  return vec(density)
end
