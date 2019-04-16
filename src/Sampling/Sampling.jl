export SamplingPattern

mutable struct SamplingPattern
  shape::Tuple
  redFac::Float64
  patParams
end

include("Simple.jl")
include("Regular.jl")
include("Vardens.jl")
include("Lines.jl")
include("PoissonDisk.jl")
include("VDPoissonDisk.jl")
include("CalibrationArea.jl")
include("PointSpreadFunction.jl")
include("Incoherence.jl")


function SamplingPattern(shape::Tuple, redFac::Float64, patFunc::AbstractString; kargs...)

if redFac < 1
  error("Reduction factor redFac must be >= 1")
end

if patFunc == "simple"
  return SamplingPattern(shape,redFac,SimplePatternParams(;kargs...))
elseif patFunc == "regular"
  return SamplingPattern(shape,redFac,RegularPatternParams(;kargs...))
elseif patFunc == "vardens"
  return SamplingPattern(shape,redFac,VardensPatternParams(;kargs...))
elseif patFunc == "lines"
  return SamplingPattern(shape,redFac,LinesPatternParams(;kargs...))
elseif patFunc == "poisson"
  return SamplingPattern(shape,redFac,PoissonDiskPatternParams(;kargs...))
elseif patFunc == "vdPoisson"
  return SamplingPattern(shape,redFac,VDPoissonDiskParams(;kargs...))
else
  error("Sample function $(patFunc) not found.")
end

end

function sample_kspace(kspace::AbstractArray,redFac::Float64,patFunc::AbstractString;kargs...)
  sample_kspace(kspace,SamplingPattern(size(kspace),redFac,patFunc;kargs...);kargs...)
end

function sample_kspace(kspace::AbstractArray,pattern::SamplingPattern;kargs...)
  patOut = sample(pattern.shape,pattern.redFac,pattern.patParams;kargs...)
  patOut = sort(patOut)
  return kspace[patOut],patOut
end

function sample_kspace(acqData::AcquisitionData,redFac::Float64,
                       patFunc::AbstractString; rand=true, profiles=false,
                       seed = 1234, kargs...)
  acqData2 = deepcopy(acqData)
  sample_kspace!(acqData2,redFac,patFunc;rand=rand,profiles=profiles,seed=seed,kargs...)
  return acqData2
end

function sample_kspace!(acqData::AcquisitionData,redFac::Float64,
                       patFunc::AbstractString; rand=true, profiles=false,
                       seed = 1234, kargs...)

  numEchoes = acqData.numEchoes
  numCoils = acqData.numCoils
  numSlices = acqData.numSlices

  idx = Vector{Array{Int64}}(undef,numEchoes)

  for echo = 1:numEchoes
    tr = trajectory(acqData,echo)
    samplingShape = ( numSamplingPerProfile(tr), numProfiles(tr) )
    pattern = SamplingPattern(samplingShape, redFac, patFunc; seed = seed, kargs...)
    patOut = sample(samplingShape,redFac,pattern.patParams; seed = seed, kargs...)
    patOut = sort(patOut)
    for slice=1:numSlices
      acqData.kdata[echo,slice] = acqData.kdata[echo,slice][patOut,:]
    end
    acqData.subsampleIndices[echo] = patOut
    rand && (seed += 1)
  end
end

function shuffle_vector(vec::Vector{T};patFunc::AbstractString="poisson",redFac::Float64=one(Float64),kargs...) where T
  shuffle_vector(vec,redFac,patFunc;kargs...)
end

function shuffle_vector(vec::Vector{T},redFac::Float64,patFunc::AbstractString;kargs...) where T
  shuffle_vector(vec,SamplingPattern(size(vec),redFac,patFunc;kargs...);kargs...)
end

function shuffle_vector(vec::Vector{T},pattern::SamplingPattern;kargs...) where T
  patOut = sample(pattern.shape,pattern.redFac,pattern.patParams;kargs...)
  return vec[patOut]
end
