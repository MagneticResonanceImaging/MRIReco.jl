export CartesianTrajectory, cartesianTrajectory2dNodes, cartesian2dDensity


"""
    CartesianTrajectory(T, numProfiles, numSamplingPerProfile
                  ; TE::T=0.0
                  , AQ::T=1.e-3
                  , kmin=(-0.5,-0.5)
                  , kmax=(0.5,0.5)
                  , kargs...)

returns a 2d cartesian trajectory.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`kmin=(-0.5,-0.5)`)            - minimum values of the covered k-space (for partial Fourier imaging)
* (`kmax=(-0.5,-0.5)`)            - maximum values of the covered k-space (for partial Fourier imaging)
"""
function CartesianTrajectory(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , kmin=(-0.5,-0.5)
                  , kmax=(0.5,0.5)
                  , kargs...) where T
  nodes = cartesian2dNodes(T, numProfiles, numSamplingPerProfile; kmin=kmin,kmax=kmax)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("Cartesian", nodes, times, T(TE), T(AQ), numProfiles, numSamplingPerProfile, 1, true, false)
end

function cartesian2dNodes(::Type{T}, numProfiles, numSamplingPerProfile
                  ; kmin=(-0.5,-0.5)
                  , kmax=(0.5,0.5)
                  , kargs...) where T

  nodes = zeros(T, 2, numSamplingPerProfile, numProfiles)
  posX = collect( -ceil(Int64, (numSamplingPerProfile-1)/2.):floor(Int64, (numSamplingPerProfile-1)/2.) ) / numSamplingPerProfile
  posY = collect( -ceil(Int64, (numProfiles-1)/2.):floor(Int64, (numProfiles-1)/2.) ) / numProfiles

  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      nodes[1,k,l] = posX[k]
      nodes[2,k,l] = posY[l]
    end
  end

  # rescale nodes according to the region covered by kmin & kmax
  nodes[1,:,:] .*= (kmax[1]-kmin[1])
  nodes[2,:,:] .*= (kmax[2]-kmin[2])
  # shift nodes to the region defined by kmin & kmax
  nodes[1,:,:] .+= 0.5*(kmin[1]+kmax[1])
  nodes[2,:,:] .+= 0.5*(kmin[2]+kmax[2])

  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end

function cartesianSubsamplingIdx(shape::NTuple{2,Int}, tr::Trajectory)
  if !isCartesian(tr)
    @error "SampledFFTOp can only be applied to Cartesian data"
  end
  
  nx,ny = shape
  numProf, numSamp = (numProfiles(tr), numSamplingPerProfile(tr))
  idxX = range(1,step=floor(Int, nx/numSamp), length=numSamp)
  idxY = range(1,step=floor(Int, ny/numProf), length=numProf)

  idx = Vector{Int}(undef,numProf*numSamp)
  linIdx1 = LinearIndices((nx,ny))
  linIdx2 = LinearIndices((numSamp,numProf))
  for j=1:numProf, i=1:numSamp
    idx[linIdx2[i,j]] = linIdx1[idxX[i],idxY[j]]
  end

  return idx
end

function isUndersampledCartTrajectory(shape::NTuple{2,Int}, tr::Trajectory)
  return shape[1]!=numSamplingPerProfile(tr) || shape[2]!=numProfiles(tr)
end

function cartesian2dDensity(::Type{T}, numProfiles::Int64, numSamplingPerProfile::Int64) where T
  density = zeros(T,numSamplingPerProfile, numProfiles) 
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      density[k,l] = 1
    end
  end
  return vec(density)
end
