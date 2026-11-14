module MRISampling

export sample, sample_kspace, sample_random, sample_regular, sample_lines,
       sample_poissondisk, sample_vdpoisson

using MRIBase
using Random

include("Regular.jl")
include("Random.jl")
include("Lines.jl")
include("PoissonDisk.jl")
include("VDPoissonDisk.jl")
include("CalibrationArea.jl")
include("PointSpreadFunction.jl")
include("Incoherence.jl")

"""
    sample(shape::NTuple{N,Int64}, redFac, patFunc::String; kargs...)

generates a `Vector{Int64}` of indices to sample an Array of of size `shape` with
a reduction factor `redFac`.

# Arguments
* `shape::NTuple{N,Int64}` - size of the Array to be sampled
* `redFac::Float64`        - subsampling factor
* `patFunc::String`        - name of the sampling function
                            ("random, "regular", "lines", "poisson" or "vdPoisson")
* `kargs...`               - additional keyword arguments of the sampling function,
                            for example `rng` and `seed`
"""
function sample(shape::NTuple{N,Int64}, redFac, patFunc::String; kargs...) where N
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

"""
    sample_kspace(data::AbstractArray, redFac::Float64, patFunc::String; kargs...)

subsamples the Array `data` with a reduction factor `redFac`
and returns both the subsampled Array (as a vector) and the sampled indices (as a vector)

# Arguments
* `data::AbstractArray`    - array to be sampled
* `redFac::Float64`        - subsampling factor
* `patFunc::String`        - name of the sampling function
                            ("random, "regular", "lines", "poisson" or "vdPoisson")
* `kargs...`               - addional keyword arguments
"""
function sample_kspace(data::AbstractArray, redFac, patFunc::AbstractString;kargs...)
  patOut = sample(size(data),redFac,patFunc;kargs...)
  patOut = sort(patOut)
  return data[patOut],patOut
end

"""
    sample_kspace(acqData::AcquisitionData,redFac::Float64,
                patFunc::AbstractString; rand=true, profiles=true,
                rng = Random.default_rng(), seed = nothing, kargs...)

subsamples the data in `acqData` with reduction factor `redFac`
and returns a new `AcquisitionData` object.

# Arguments
* `acqData::AcquisitionDatay` - AcquisitionData to be sampled
* `redFac::Float64`           - subsampling factor
* `patFunc::String`           - name of the sampling function
                                ("random, "regular", "lines", "poisson" or "vdPoisson")
* (`rand=true`)               - use different patterns for the different contrasts
* (`profiles=true`)           - sample complete profiles
* (`rng=Random.default_rng()`) - random number generator
* (`seed=nothing`)            - seed for `rng`; with `nothing`, `rng` is used in its current state
* `kargs...`                  - addional keyword arguments
"""
function sample_kspace(acqData::AcquisitionData, redFac,
                       patFunc::AbstractString; rand=true, profiles=true,
                       rng = Random.default_rng(), seed = nothing, kargs...)
  acqData2 = deepcopy(acqData)
  sample_kspace!(acqData2,redFac,patFunc;rand=rand,profiles=profiles,rng=rng,seed=seed,kargs...)
  return acqData2
end

"""
    sample_kspace(acqData::AcquisitionData,redFac::Float64,
                patFunc::AbstractString; rand=true, profiles=true,
                rng = Random.default_rng(), seed = nothing, kargs...)

subsamples the data in acqData object with reduction factor `redFac`. Performs the same
subsampling as `sample_kspace` but operates in-place on acqData.
"""
function sample_kspace!(acqData::AcquisitionData, redFac,
                       patFunc::AbstractString; rand=true, profiles=true,
                       rng = Random.default_rng(), seed = nothing, kargs...)

  numContr, numSl, numRep = numContrasts(acqData), numSlices(acqData), numRepetitions(acqData)

  for echo = 1:numContr
    tr = trajectory(acqData,echo)
    if ndims(tr)==2
      samplingShape = ( numSamplingPerProfile(tr), numProfiles(tr) )
    else
      samplingShape = ( numSamplingPerProfile(tr), numProfiles(tr), numSlices(tr) )
    end
    # without a seed, the same pattern for all contrasts comes from the same state of the generator,
    # which only the last contrast advances
    patRng = (rand || !isnothing(seed) || echo == numContr) ? rng : copy(rng)
    # only sample full profiles
    if profiles
      patOut = sample(samplingShape,redFac,"lines"; sampleFunc=patFunc, rng=patRng, seed=seed, kargs...)
    else
      patOut = sample(samplingShape,redFac,patFunc; rng=patRng, seed=seed, kargs...)
    end
    patOut = sort(patOut)
    for slice=1:numSl, rep=1:numRep
      acqData.kdata[echo,slice, rep] = acqData.kdata[echo,slice,rep][patOut,:]
    end
    acqData.subsampleIndices[echo] = patOut
    rand && !isnothing(seed) && (seed += 1)
  end
end


end # module
