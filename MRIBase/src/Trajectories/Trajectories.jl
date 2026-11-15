export Trajectory, trajectory,
       kspaceNodes, echoTime, acqTimePerProfile, readoutTimes,
       numSamplingPerProfile, numProfiles, numSlices, samplingFreq,
       scale, isCartesian, isCircular

import Base: string

"""
struct describing a trajectory with element type `T`

# Fields
* `name::String`                  - name of the trajectory
* `nodes::Matrix{T}`              - sampling locations in k-space.
                                    (1.dim <-> dimensions of k-space, 2. dim <-> sampling points)
* `times::Vector{T}`              - sampling times in s
* `TE::T`                         - echo time in s
* `AQ::T`                         - readout duration in s (per profile)
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* `numSlices::Int64`              - number of slices (for 3d trajectories)
* `cartesian::Bool`               - true if sampling points lie on a cartesian grid
* `circular::Bool`                - true if kspace is covered in a circular domain
"""
mutable struct Trajectory{T}
  name::String
  nodes::Matrix{T}
  times::Vector{T}
  TE::T
  AQ::T
  numProfiles::Int64
  numSamplingPerProfile::Int64
  numSlices::Int64
  cartesian::Bool
  circular::Bool
end

"""
    Trajectory(nodes::AbstractMatrix{T}, numProfiles::Int64, numSamplingPerProfile;
               times=nothing, TE=0.0, AQ=1.e-3, numSlices::Int64=1,
               cartesian::Bool=false, circular::Bool=false) where T <: AbstractFloat

creates a trajectory named `"Custom"` from the k-space positions `nodes`, a `Matrix{T}` with one row per dimension
and one column per sample, ordered profile by profile and, in a stacked 3d trajectory, partition by partition.
Without `times`, the samples of each profile are spread over the readout duration `AQ`, starting at `TE`.

# Arguments
* `nodes::Matrix{T}`          - k-space positions between -0.5 and 0.5 (other matrix types, such as a transpose,
                                have to be converted with `Matrix` first)
* `numProfiles::Int64`        - number of profiles per partition
* `numSamplingPerProfile`     - number of sampling points per profile
* (`times=nothing`)           - sampling times in s, a `Vector{T}` with one time per node
* (`TE=0.0`)                  - echo time in s
* (`AQ=1.e-3`)                - readout duration in s (per profile)
* (`numSlices::Int64=1`)      - number of partitions (slices) of a stacked 3d trajectory
* (`cartesian::Bool=false`)   - true if the nodes lie on the cartesian grid in the order of `CartesianTrajectory`
                                (`CartesianTrajectory3D` in 3d): without a correction map, the reconstruction and
                                the default simulation then use FFTs, which ignore the nodes
* (`circular::Bool=false`)    - true if the nodes cover a circular region of k-space, as in radial and spiral
                                trajectories: the reconstructions `"direct"`, `"standard"` and `"multiCoil"` then set
                                the image outside a circle (a sphere in 3d) whose diameter is the largest side of the
                                image to zero, and with a correction map, the encoding operator ignores the image
                                outside it
"""
function Trajectory(nodes::AbstractMatrix{T}, numProfiles::Int64, numSamplingPerProfile;
                    times=nothing, TE=0.0, AQ=1.e-3, numSlices::Int64=1,
                    cartesian::Bool=false, circular::Bool=false) where T <: AbstractFloat

  TE_ = T(TE)
  AQ_ = T(AQ)

  if times != nothing
    ttimes = times
  else
    ttimes = readoutTimes(T,numProfiles,numSamplingPerProfile, numSlices; TE=TE_, AQ=AQ_)
  end

  return Trajectory("Custom", nodes, ttimes, TE_, AQ_, numProfiles, numSamplingPerProfile, numSlices, cartesian, circular)
end

Base.vec(tr::Trajectory) = [tr]
Base.vec(tr::Vector{Trajectory}) = tr
Base.size(tr::Trajectory) = size(kspaceNodes(tr))
Base.size(tr::Trajectory,i::Int) = size(kspaceNodes(tr),i)

"""
    scale(tr::Trajectory, factor::Real)

returns a copy of `tr` with the nodes divided by `factor`
"""
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
    trajectory(::Type{T}, trajName::AbstractString, numProfiles::Int, numSamplingPerProfile::Int;
               numSlices::Int64=1, TE::AbstractFloat=0.0, AQ::AbstractFloat=1.e-3, kargs...) where T

is a factory method to construct a trajectory with element type `T`, for example `Float64`, from its name `trajName`

# Arguments
* `T`                             - element type of the nodes and times
* `trajName::AbstractString`      - name of the trajectory
* `numProfiles::Int`              - number of profiles
* `numSamplingPerProfile::Int`    - number of sampling points per profile
* (`numSlices::Int64=1`)          - number of slices (for 3d trajectories)
* (`TE::AbstractFloat=0.0`)       - echo time in s
* (`AQ::AbstractFloat=1.e-3`)     - readout duration in s (per profile)
* `kargs...`                      - additional keyword arguments, passed on to the constructor of the trajectory
"""
function trajectory(::Type{T}, trajName::AbstractString, numProfiles::Int, numSamplingPerProfile::Int; numSlices::Int64=1, 
         TE::AbstractFloat=0.0, AQ::AbstractFloat=1.e-3, kargs...) where T

  TE_ = T(TE)
  AQ_ = T(AQ)

  if trajName == "Spiral"
    tr = SpiralTrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
  elseif trajName == "Radial"
    tr = RadialTrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
  elseif trajName == "Cartesian"
    tr = CartesianTrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
  elseif trajName == "EPI"
    tr = EPITrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
  elseif trajName == "OneLine"
    tr = OneLine2dTrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
  elseif trajName == "SpiralVarDens"
    tr = SpiralTrajectoryVarDens(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
elseif trajName == "SpiralPerturbed"
    tr = SpiralPerturbedTrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
  elseif trajName == "SpiralDualDens"
    tr = SpiralTrajectoryDualDens(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
  elseif trajName == "Cartesian3D"
    tr = CartesianTrajectory3D(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, numSlices=numSlices, kargs...)
  elseif trajName == "StackOfStars"
    tr = StackOfStarsTrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, numSlices=numSlices, kargs...)
  elseif trajName == "Kooshball"
    tr = KooshballTrajectory(T, numProfiles, numSamplingPerProfile; TE=TE_, AQ=AQ_, kargs...)
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

""" `numProfiles(tr::Trajectory)` returns the number of profiles of a trajectory """
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
""" `ndims(tr::Trajectory)` returns the number of dimensions of a trajectory """
Base.ndims(tr::Trajectory) = size(tr.nodes,1)

""" `kspaceNodes(tr::Trajectory)` returns the kspace sampling points of a trajectory """
kspaceNodes(tr::Trajectory) = tr.nodes

""" `readoutTimes(tr::Trajectory)` returns the readoutTimes for the sampling points of a trajectory """
readoutTimes(tr::Trajectory) = tr.times

""" `profileNodes(tr::Trajectory,prof::Int,slice::Int)` returns kspace nodes for a given profile """
function profileNodes(tr::Trajectory,prof::Int,slice::Int)
  nodes = reshape(tr.nodes,:,tr.numSamplingPerProfile,tr.numProfiles, tr.numSlices)
  return nodes[:,:,prof,slice]
end

function readoutTimes(::Type{T}, numProfiles, numSamplingPerProfile, numSlices=1; TE=0.0, AQ=1.e-3) where T
  times = zeros(T, numSamplingPerProfile, numProfiles, numSlices)
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
