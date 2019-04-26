export profileIdx, sample

include("Regular.jl")
include("Random.jl")
include("Lines.jl")
include("PoissonDisk.jl")
include("VDPoissonDisk.jl")
include("CalibrationArea.jl")
include("PointSpreadFunction.jl")
include("Incoherence.jl")

function sample(shape::NTuple{N,Int64}, redFac::Float64, patFunc::String; kargs...) where N
  if redFac < 1
    error("Reduction factor redFac must be >= 1")
  end

  if patFunc == "random"
    return sample_random(shape,redFac;kargs...)
  elseif patFunc == "regular"
    return sample_regular(shape,redFac;kargs...)
  elseif patFunc == "lines"
    return sample_lines(shape,redFac;kargs...)
  elseif patFunc == "poisson"
    return sample_poissondisk(shape,redFac;kargs...)
  elseif patFunc == "vdPoisson"
    return sample_vdpoisson(shape,redFac;kargs...)
  else
    @error "Sample function $(patFunc) not found."
  end
end

function sample_kspace(data::AbstractArray,redFac::Float64,patFunc::AbstractString;kargs...)
  patOut = sample(size(data),redFac,patFunc;kargs...)
  patOut = sort(patOut)
  return data[patOut],patOut
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
  numSl = acqData.numSlices

  idx = Vector{Array{Int64}}(undef,numEchoes)

  for echo = 1:numEchoes
    tr = trajectory(acqData,echo)
    if dims(tr)==2
      samplingShape = ( numSamplingPerProfile(tr), numProfiles(tr) )
    else
      samplingShape = ( numSamplingPerProfile(tr), numProfiles(tr), numSlices(tr) )
    end
    # only sample full profiles
    patOut = sample(samplingShape,redFac,"lines"; sampleFunc=patFunc, seed=seed, kargs...)
    patOut = sort(patOut)
    for slice=1:numSl
      acqData.kdata[echo,slice] = acqData.kdata[echo,slice][patOut,:]
    end
    acqData.subsampleIndices[echo] = patOut
    rand && (seed += 1)
  end
end

function profileIdx(acqData::AcquisitionData,contr::Int)
  tr = trajectory(acqData,contr)
  numSamp = numSamplingPerProfile(tr)
  numProf = div(length(acqData.subsampleIndices[contr]),numSamp)
  idx = [div(acqData.subsampleIndices[contr][numSamp*(prof-1)+1],numSamp)+1 for prof=1:numProf]
end
