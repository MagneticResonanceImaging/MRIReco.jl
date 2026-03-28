export StandardReconstruction
@reconstruction mutable struct StandardReconstruction{P<:AbstractStandardParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end

export MultiEchoReconstruction  
@reconstruction mutable struct MultiEchoReconstruction{P<:AbstractMultiEchoParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end

export MultiCoilReconstruction
@reconstruction mutable struct MultiCoilReconstruction{P<:AbstractMultiCoilParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end

export MultiCoilMultiEchoReconstruction
@reconstruction mutable struct MultiCoilMultiEchoReconstruction{P<:AbstractMultiCoilMultiEchoParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end

export SubspaceReconstruction
@reconstruction mutable struct SubspaceReconstruction{P<:AbstractSubspaceParameters} <: AbstractIterativeMRIRecoAlgorithm
    @parameter parameter::P
end

#(params::DirectMRIRecoParameter{Nothing})(algo::AbstractDirectMRIRecoAlgorithm, acqData::AcquisitionData) = params(algo, acqData, encodingSize(acqData))
#(params::DirectMRIRecoParameter)(algo::AbstractDirectMRIRecoAlgorithm, acqData::AcquisitionData) = params(algo, acqData, params.reconSize)
