export MultiCoilReconstruction

export AbstractMultiCoilParameters
abstract type AbstractMultiCoilParameters <: AbstractIterativeRecoParameters end

@reconstruction mutable struct MultiCoilReconstruction{P<:AbstractMultiCoilParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end