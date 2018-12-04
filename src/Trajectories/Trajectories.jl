export AbstractTrajectory,
       Abstract2DTrajectory, Abstract3DTrajectory, CartesianTrajectory,
       StackedTrajectory, Real3DTrajectory,
       kspaceNodes, string, echoTime, acqTimePerProfile, readoutTimes,
       numSamplingPerProfile, numProfiles, findCenters

import Base: string

abstract type AbstractTrajectory end

abstract type Abstract2DTrajectory <: AbstractTrajectory end
abstract type CartesianTrajectory <: Abstract2DTrajectory end
abstract type Abstract3DTrajectory <: AbstractTrajectory end

abstract type StackedTrajectory <: Abstract3DTrajectory end
abstract type Real3DTrajectory <: Abstract3DTrajectory end

include("2D/Cartesian2D.jl")
include("2D/Radial2D.jl")
include("2D/Spiral2D.jl")
include("2D/OneLine2D.jl")
include("2D/Spiral2DVariableDens.jl")
include("2D/EPI.jl")

include("3D/Kooshball.jl")
include("3D/Cartesian3D.jl")
include("3D/StackofStars.jl")

include("UndersampledTrajectory.jl")
include("CustomTrajectory.jl")

#### Factory method ###

# This dispatches on the file extension and automatically
# generates the correct type
function trajectory(trajName::AbstractString, numProfiles::Int, numSamplingPerProfile::Int; numSlices::Int64=1, nodes=nothing, kargs...)
  if trajName == "Spiral"
    return SpiralTrajectory(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "Radial"
    return RadialTrajectory(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "Cartesian"
    return SimpleCartesianTrajectory(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "EPI"
    return EPITrajectory(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "OneLine"
    return OneLine2DTrajectory(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "SpiralVarDens"
    return SpiralTrajectoryVarDens(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "Cartesian3D"
    return CartesianTrajectory3D(numProfiles, numSamplingPerProfile; numSlices=numSlices, kargs...)
  elseif trajName == "StackOfStars"
    return StackOfStarsTrajectory(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "Kooshball"
    return KooshballTrajectory(numProfiles, numSamplingPerProfile; kargs...)
  elseif trajName == "Custom2D"
    return CustomTrajectory(numProfiles, numSamplingPerProfile, nodes; kargs...)
  elseif trajName == "Custom3D"
    return CustomTrajectory(numProfiles, numSamplingPerProfile, numSlices, nodes; kargs...)
  else
    error("The trajectory $trajName is not yet supported!")
  end
  # make this function type stable...
  return CartesianTrajectory(1,1)
end


## some general methods that might be overloaded by the real implementations

echoTime(tr::AbstractTrajectory) = tr.TE
acqTimePerProfile(tr::AbstractTrajectory) = tr.AQ
numProfiles(tr::AbstractTrajectory) = tr.numProfiles
numSamplingPerProfile(tr::AbstractTrajectory) = tr.numSamplingPerProfile

numSlices(tr::AbstractTrajectory) = tr.numSlices
isCircular(tr::AbstractTrajectory) = false

# this function calculates the readout time for each sampling point on the
# kspace tracectory
function readoutTimes(tr::Abstract2DTrajectory)
  times = zeros(numSamplingPerProfile(tr), numProfiles(tr),)
  for l = 1:numProfiles(tr)
    for k = 1:numSamplingPerProfile(tr)
      times[k,l] = (echoTime(tr) +
                   acqTimePerProfile(tr)*(k-1)/numSamplingPerProfile(tr) )
    end
  end
  return times
end

# this function calculates the readout time for each sampling point on the
# kspace tracectory
function readoutTimes(tr::StackedTrajectory)
  times = zeros(numSamplingPerProfile(tr), numProfiles(tr), numSlices(tr))
  for j = 1:numSlices(tr)
    for l = 1:numProfiles(tr)
      for k = 1:numSamplingPerProfile(tr)
        times[k,l,j] = (echoTime(tr) +
                   acqTimePerProfile(tr)*(k-1)/numSamplingPerProfile(tr) )
      end
    end
  end
  return times
end

function readoutTimes(tr::Real3DTrajectory)
  times = zeros(numSamplingPerProfile(tr), numProfiles(tr),)
  for l = 1:numProfiles(tr)
    for k = 1:numSamplingPerProfile(tr)
      times[k,l] = (echoTime(tr) +
                   acqTimePerProfile(tr)*(k-1)/numSamplingPerProfile(tr) )
    end
  end
  return times
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

function findCenters(tr::Abstract2DTrajectory)
  knodes = kspaceNodes(tr)
  temp = find(x->x==0,knodes[1,:])

  zeroInd=Int64[]

  for i in temp

    if knodes[2,i] == 0
      push!(zeroInd,i)
    end

  end

  return zeroInd

end

function findCenters(tr::Abstract3DTrajectory)
  knodes = kspaceNodes(tr)
  temp = find(x->x==0,knodes[1,:])

  zeroIndY=Int64[]
  zeroInd = Int64[]

  for i in temp
    if knodes[2,i] == 0
      push!(zeroIndY,i)
    end
  end

  for i in zeroIndY
    if knodes[3,i] == 0
      push!(zeroInd,i)
    end
  end

  return zeroInd

end

function findCenters(tr::StackOfStarsTrajectory)
  knodes = kspaceNodes(tr)
  temp = find(x->x==0,knodes[1,:])

  zeroInd = Int64[]

  for i in temp
    if knodes[2,i] == 0
      push!(zeroInd,i)
    end
  end

  return zeroInd

end

function findCenterNeighbours(tr::CartesianTrajectory,range::Int64)
  knodes = kspaceNodes(tr)
  pRange = minimum(knodes[1,find(x->x>0,knodes[1,:])])
  mRange = maximum(knodes[1,find(x->x<0,knodes[1,:])])

  idx=Int64[]

  zeroInd = findCenters(tr)

  push!(idx,zeroInd[1])


  for j=1:range

    temp = find(x->x==0,knodes[1,:])

    for i in temp
      for k=range*mRange:pRange:range*pRange
        if k!=0
          if knodes[2,i] == k
            push!(idx,i)
          end
        end
      end
    end

    temp = find(x->x==j*pRange,knodes[1,:])

    for i in temp
      for k=range*mRange:pRange:range*pRange
        if knodes[2,i] == k
          push!(idx,i)
        end
      end
    end

    temp = find(x->x==j*mRange,knodes[1,:])

    for i in temp
      for k=range*mRange:pRange:range*pRange
        if knodes[2,i] == k
          push!(idx,i)
        end
      end
    end

  end

  return idx

end

function improvedTrajPlot(tr::Abstract2DTrajectory)
  numProfiles = tr.numProfiles
  numSamplePerProfile = tr.numSamplingPerProfile
  nodes = kspaceNodes(tr)
  figure()
  for i=1:numProfiles
     indices = 1 + (i-1)*numSamplePerProfile : i*numSamplePerProfile
     plot(nodes[1,indices],nodes[2,indices],color="blue" )
     plot(nodes[1,indices],nodes[2,indices],".",color="red", )
  end
end
