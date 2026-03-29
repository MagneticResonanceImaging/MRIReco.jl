include("SolverParameters.jl")
include("WeightingParameters.jl")

export AbstractIterativeRecoParameters
abstract type AbstractIterativeRecoParameters <: AbstractMRIRecoParameters end

export AbstractStandardParameters
abstract type AbstractStandardParameters <: AbstractIterativeRecoParameters end
export AbstractMultiEchoParameters  
abstract type AbstractMultiEchoParameters <: AbstractIterativeRecoParameters end
export AbstractMultiCoilParameters
abstract type AbstractMultiCoilParameters <: AbstractIterativeRecoParameters end
export AbstractMultiCoilMultiEchoParameters
abstract type AbstractMultiCoilMultiEchoParameters <: AbstractIterativeRecoParameters end
export AbstractSubspaceParameters
abstract type AbstractSubspaceParameters <: AbstractIterativeRecoParameters end

#=
struct StandardParameters <: AbstractStandardParameters
    encoding::EncodingParameters
    solver::SolverParameter
end
# src/Algorithms/Parameters/MultiEchoParameters.jl  
struct MultiEchoParameters <: AbstractMultiEchoParameters
    encoding::EncodingParameters
    solver::SolverParameter
    numEchoes::Int
end
# src/Algorithms/Parameters/MultiCoilParameters.jl
struct MultiCoilParameters <: AbstractMultiCoilParameters
    encoding::EncodingParameters
    solver::SolverParameter
    senseMaps::AbstractArray
    L_inv::Union{LowerTriangular, Nothing}
end
# src/Algorithms/Parameters/MultiCoilMultiEchoParameters.jl
struct MultiCoilMultiEchoParameters <: AbstractMultiCoilMultiEchoParameters
    encoding::EncodingParameters
    solver::SolverParameter
    senseMaps::AbstractArray
    L_inv::Union{LowerTriangular, Nothing}
    numEchoes::Int
end
# src/Algorithms/Parameters/SubspaceParameters.jl
struct SubspaceParameters <: AbstractSubspaceParameters
    encoding::EncodingParameters
    solver::SolverParameter
    basis::AbstractMatrix
    numBasis::Int
    senseMaps::AbstractArray
    L_inv::Union{LowerTriangular, Nothing}
end
=#