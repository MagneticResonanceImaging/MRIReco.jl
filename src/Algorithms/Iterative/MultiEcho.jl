export MultiEchoReconstruction
export AbstractMultiEchoParameters  

abstract type AbstractMultiEchoParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct MultiEchoReconstruction{P<:AbstractMultiEchoParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end