export SubspaceReconstruction

export AbstractSubspaceParameters
abstract type AbstractSubspaceParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct SubspaceReconstruction{P<:AbstractSubspaceParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end