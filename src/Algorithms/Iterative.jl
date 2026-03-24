abstract type AbstractIterativeRecoParameters <: AbstractMRIRecoParameters end

@reconstruction mutable struct IterativeMRIReconstruction{P <: AbstractIterativeRecoParameters} <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::P
end
