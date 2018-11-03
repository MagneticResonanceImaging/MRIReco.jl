export UndersampledTrajectory, readoutTimes

mutable struct UndersampledTrajectory <: AbstractTrajectory
  trajectory::AbstractTrajectory
  subsampleIndices::Vector{Int64}
end

#
# return only the nodes of  trajectory which are sampled  (specified by tr.subsampleIndices)
#
function kspaceNodes(tr::UndersampledTrajectory)
  nodes = kspaceNodes(tr.trajectory)
  return nodes[:,tr.subsampleIndices]
end

function kspaceDensity(tr::UndersampledTrajectory)
  error("Not implemented!")
end

readoutTimes(tr::UndersampledTrajectory) = readoutTimes(tr.trajectory)[tr.subsampleIndices]
echoTime(tr::UndersampledTrajectory) = echoTime(tr.trajectory)
