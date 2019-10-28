export Trajectory,trajectory,
       kspaceNodes, echoTime, acqTimePerProfile, readoutTimes,
       numSamplingPerProfile, numProfiles, numSlices, findCenters, samplingFreq,
       scale

import Base: string

"""
struct describing a trajectory

# Fields
* `name::String`                  - name of the trajectory
* `nodes::Matrix{Float64}`        - sampling locations in k-space.
                                    (1.dim <-> dimensions of k-space, 2. dim <-> sampling points)
* `times::Vector{Float64}`        - sampling times in s
* `TE::Float64`                   - echo time in s
* `AQ::Float64``                  - readout duration in s (per profile)
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* `numSlices::Int64`              - number of slices (for 3d trajectories)
* `cartesian::Bool`               - true if sampling points lie on a cartesian grid
* `circular::Bool`                - true if kspace is covered in a circular domain
"""
mutable struct Trajectory
  name::String
  nodes::Matrix{Float64}
  times::Vector{Float64}
  TE::Float64
  AQ::Float64
  numProfiles::Int64
  numSamplingPerProfile::Int64
  numSlices::Int64
  cartesian::Bool
  circular::Bool
end

function Trajectory(nodes::AbstractMatrix{T}, numProfiles::Int64, numSamplingPerProfile;
                    times=nothing, TE::Float64=0.0, AQ::Float64=1.e-3, numSlices::Int64=1,
                    cartesian::Bool=false, circular::Bool=false) where T <: AbstractFloat
  if times != nothing
    ttimes = readoutTimes
  else
    ttimes = readoutTimes(numProfiles,numSamplingPerProfile; TE=TE, AQ=AQ)
  end

  return Trajectory("Custom", Float64.(nodes), ttimes, TE, AQ, numProfiles, numSamplingPerProfile, numSlices, cartesian, circular)
end

Base.vec(tr::Trajectory) = [tr]
Base.vec(tr::Vector{Trajectory}) = tr
Base.size(tr::Trajectory) = size(kspaceNodes(tr))
Base.size(tr::Trajectory,i::Int) = size(kspaceNodes(tr),i)
function scale(tr::Trajectory, factor::Real)
  trC = deepcopy(tr)
  trC.nodes ./= factor
  return trC
end

include("2D/Cartesian2D.jl")
include("2D/Radial2D.jl")
include("2D/Spiral2D.jl")
include("2D/OneLine2D.jl")
include("2D/Spiral2DVariableDens.jl")
include("2D/Spiral2DDualDens.jl")
include("2D/Spiral2DPerturbed.jl")
include("2D/EPI.jl")

include("3D/Kooshball.jl")
include("3D/Cartesian3D.jl")
include("3D/StackofStars.jl")

#### Factory method ###

# This dispatches on the file extension and automatically
# generates the correct type
"""
    trajectory(trajName::AbstractString, numProfiles::Int, numSamplingPerProfile::Int; numSlices::Int64=1, TE::Float64=0.0, AQ::Float64=1.e-3, kargs...)

is a factory method to construct a trajectory from its `name`

# Arguments
* `name::String`                  - name of the trajectory
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`numSlices::Int64=1`)          - number of slices (for 3d trajectories)
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* `kargs...`                      - addional keyword arguments
"""
function trajectory(trajName::AbstractString, numProfiles::Int, numSamplingPerProfile::Int; numSlices::Int64=1, TE::Float64=0.0, AQ::Float64=1.e-3, kargs...)
  if trajName == "Spiral"
    tr = SpiralTrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  elseif trajName == "Radial"
    tr = RadialTrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  elseif trajName == "Cartesian"
    tr = CartesianTrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  elseif trajName == "EPI"
    tr = EPITrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  elseif trajName == "OneLine"
    tr = OneLine2dTrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  elseif trajName == "SpiralVarDens"
    tr = SpiralTrajectoryVarDens(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
elseif trajName == "SpiralPerturbed"
    tr = SpiralPerturbedTrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  elseif trajName == "SpiralDualDens"
    tr = SpiralTrajectoryDualDens(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  elseif trajName == "Cartesian3D"
    tr = CartesianTrajectory3D(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, numSlices=numSlices, kargs...)
  elseif trajName == "StackOfStars"
    tr = StackOfStarsTrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, numSlices=numSlices, kargs...)
  elseif trajName == "Kooshball"
    tr = KooshballTrajectory(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ, kargs...)
  else
    @error "The trajectory $trajName is not yet supported!"
  end

  # return Trajectory(trajName, nodes, TE, AQ, numProfiles, numSamplingPerProfile, numSlices, false, true)
  return tr
end


## some general methods that might be overloaded by the real implementations
""" `string(tr::Trajectory)` returns the name of a trajectory """
string(tr::Trajectory) = tr.name

""" `echoTime(tr::Trajectory)` returns the echo time of a trajectory """
echoTime(tr::Trajectory) = tr.TE

""" `acqTimePerProfile(tr::Trajectory)` returns the acquisition time per profile of a trajectory """
acqTimePerProfile(tr::Trajectory) = tr.AQ

""" `string(tr::Trajectory)` returns the name of a trajectory """
numProfiles(tr::Trajectory) = tr.numProfiles

""" `numSamplingPerProfile(tr::Trajectory)` returns the number of samples per profile of a trajectory """
numSamplingPerProfile(tr::Trajectory) = tr.numSamplingPerProfile

""" `numSlices(tr::Trajectory)` returns the number of slices of a trajectory """
numSlices(tr::Trajectory) = tr.numSlices

""" `isCircular(tr::Trajectory)` returns whether the trajectory has circular k-space coverage """
isCircular(tr::Trajectory) = tr.circular

""" `isCartesian(tr::Trajectory)` returns whether the trajectory nodes lie on a cartesian grid"""
isCartesian(tr::Trajectory) = tr.cartesian

# 2d or 3d encoding
""" `dims(tr::Trajectory)` returns the number of dimensions of a trajectory """
dims(tr::Trajectory) = size(tr.nodes,1)

""" `kspaceNodes(tr::Trajectory)` returns the kspace sampling points of a trajectory """
kspaceNodes(tr::Trajectory) = tr.nodes

""" `readoutTimes(tr::Trajectory)` returns the readoutTimes for the sampling points of a trajectory """
readoutTimes(tr::Trajectory) = tr.times

""" `profileNodes(tr::Trajectory,prof::Int,slice::Int)` returns kspace nodes for a given profile """
function profileNodes(tr::Trajectory,prof::Int,slice::Int)
  nodes = reshape(tr.nodes,:,tr.numSamplingPerProfile,tr.numProfiles, tr.numSlices)
  return nodes[:,:,prof,slice]
end

function readoutTimes(numProfiles, numSamplingPerProfile, numSlices=1; TE=0.0, AQ=1.e-3)
  times = zeros(numSamplingPerProfile, numProfiles, numSlices)
  for j = 1:numSlices
    for l = 1:numProfiles
      for k = 1:numSamplingPerProfile
        times[k,l,j] = TE + AQ*(k-1)/numSamplingPerProfile
      end
    end
  end
  return vec(times)
end


function zax(n::Int64,N::Int64)
  return (2 .*n - N -1.)/N
end

function xax(n::Int64,N::Int64)
  return cos(sqrt(N*pi)*asin(zax(n,N)))*sqrt(1.0-zax(n,N)^2.)
end

function yax(n::Int64,N::Int64)
  return sin(sqrt(N*pi)*asin(zax(n,N)))*sqrt(1.0-zax(n,N)^2.)
end

# function findCenters(tr::Abstract2DTrajectory)
#   knodes = kspaceNodes(tr)
#   temp = find(x->x==0,knodes[1,:])
#
#   zeroInd=Int64[]
#
#   for i in temp
#
#     if knodes[2,i] == 0
#       push!(zeroInd,i)
#     end
#
#   end
#
#   return zeroInd
#
# end
#
# function findCenters(tr::Trajectory)
#   knodes = kspaceNodes(tr)
#   temp = find(x->x==0,knodes[1,:])
#
#   zeroIndY=Int64[]
#   zeroInd = Int64[]
#
#   for i in temp
#     if knodes[2,i] == 0
#       push!(zeroIndY,i)
#     end
#   end
#
#   for i in zeroIndY
#     if knodes[3,i] == 0
#       push!(zeroInd,i)
#     end
#   end
#
#   return zeroInd
#
# end
#
# function findCenters(tr::StackOfStarsTrajectory)
#   knodes = kspaceNodes(tr)
#   temp = find(x->x==0,knodes[1,:])
#
#   zeroInd = Int64[]
#
#   for i in temp
#     if knodes[2,i] == 0
#       push!(zeroInd,i)
#     end
#   end
#
#   return zeroInd
#
# end
#
# function findCenterNeighbours(tr::CartesianTrajectory,range::Int64)
#   knodes = kspaceNodes(tr)
#   pRange = minimum(knodes[1,find(x->x>0,knodes[1,:])])
#   mRange = maximum(knodes[1,find(x->x<0,knodes[1,:])])
#
#   idx=Int64[]
#
#   zeroInd = findCenters(tr)
#
#   push!(idx,zeroInd[1])
#
#
#   for j=1:range
#
#     temp = find(x->x==0,knodes[1,:])
#
#     for i in temp
#       for k=range*mRange:pRange:range*pRange
#         if k!=0
#           if knodes[2,i] == k
#             push!(idx,i)
#           end
#         end
#       end
#     end
#
#     temp = find(x->x==j*pRange,knodes[1,:])
#
#     for i in temp
#       for k=range*mRange:pRange:range*pRange
#         if knodes[2,i] == k
#           push!(idx,i)
#         end
#       end
#     end
#
#     temp = find(x->x==j*mRange,knodes[1,:])
#
#     for i in temp
#       for k=range*mRange:pRange:range*pRange
#         if knodes[2,i] == k
#           push!(idx,i)
#         end
#       end
#     end
#
#   end
#
#   return idx
#
# end

function improvedTrajPlot(tr::Trajectory)
  numProfiles = numProfiles(tr)
  numSamplePerProfile = numSamplingPerProfile(tr)
  nodes = kspaceNodes(tr)
  figure()
  for i=1:numProfiles
     indices = 1 + (i-1)*numSamplePerProfile : i*numSamplePerProfile
     plot(nodes[1,indices],nodes[2,indices],color="blue" )
     plot(nodes[1,indices],nodes[2,indices],".",color="red", )
  end
end
