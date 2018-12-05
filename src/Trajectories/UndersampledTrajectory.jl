export UndersampledTrajectory, readoutTimes

mutable struct UndersampledTrajectory2D <: Abstract2DTrajectory
  trajectory::Abstract2DTrajectory
  subsampleIndices::Vector{Int64}
end

mutable struct UndersampledTrajectory3D <: Abstract3DTrajectory
  trajectory::Abstract3DTrajectory
  subsampleIndices::Vector{Int64}
end

UndersampledTrajectory(tr::Abstract2DTrajectory, subsampleIndices::Vector{Int64}) = UndersampledTrajectory2D(tr, subsampleIndices)
UndersampledTrajectory(tr::Abstract3DTrajectory, subsampleIndices::Vector{Int64}) = UndersampledTrajectory3D(tr, subsampleIndices)

#
# return only the nodes of  trajectory which are sampled  (specified by tr.subsampleIndices)
#
function kspaceNodes(tr::T) where T <: Union{UndersampledTrajectory2D,UndersampledTrajectory3D}
  nodes = kspaceNodes(tr.trajectory)
  return nodes[:,tr.subsampleIndices]
end

function kspaceDensity(tr::T) where T <: Union{UndersampledTrajectory2D,UndersampledTrajectory3D}
  error("Not implemented!")
end

readoutTimes(tr::T) where T <: Union{UndersampledTrajectory2D,UndersampledTrajectory3D} = readoutTimes(tr.trajectory)[tr.subsampleIndices]
echoTime(tr::T) where T <: Union{UndersampledTrajectory2D,UndersampledTrajectory3D} = echoTime(tr.trajectory)
isCircular(tr::T) where T <: Union{UndersampledTrajectory2D,UndersampledTrajectory3D} = isCircular(tr.trajectory)
