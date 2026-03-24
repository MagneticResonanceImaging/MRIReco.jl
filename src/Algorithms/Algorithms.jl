export AbstractMRIRecoAlgorithm
abstract type AbstractMRIRecoAlgorithm <: AbstractImageReconstructionAlgorithm end

export AbstractDirectMRIRecoAlgorithm
abstract type AbstractDirectMRIRecoAlgorithm <: AbstractMRIRecoAlgorithm end

export AbstractIterativeMRIRecoAlgorithm
abstract type AbstractIterativeMRIRecoAlgorithm <: AbstractMRIRecoAlgorithm end

export AbstractMRIRecoParameters
abstract type AbstractMRIRecoParameters <: AbstractImageReconstructionParameters end

export MRIRecoStyle
struct MRIRecoStyle <: CustomPlanStyle end

include("Storage.jl")
include("Direct.jl")
include("Iterative.jl")
include("Parameters/Parameters.jl")