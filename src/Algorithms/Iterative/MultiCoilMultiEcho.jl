export MultiCoilMultiEchoReconstruction

export AbstractMultiCoilMultiEchoParameters
abstract type AbstractMultiCoilMultiEchoParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct MultiCoilMultiEchoReconstruction{P<:AbstractMultiCoilMultiEchoParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end