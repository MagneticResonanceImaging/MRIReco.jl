export UndersampledTrajectory, readoutTimes

mutable struct UndersampledTrajectory <: AbstractTrajectory
  trajectory::AbstractTrajectory
  idx::Vector{Int64}
end

#
# return only the nodes of  trajectory which are sampled  (specified by tr.idx)
#
function kspaceNodes(tr::UndersampledTrajectory)
  nodes = kspaceNodes(tr.trajectory)
  return nodes[:,tr.idx]
end

function kspaceDensity(tr::UndersampledTrajectory)
  error("Not implemented!")
end

readoutTimes(tr::UndersampledTrajectory) = readoutTimes(tr.trajectory)[tr.idx]
echoTime(tr::UndersampledTrajectory) = echoTime(tr.trajectory)
