"""
    sample_vdpoisson(shape::Tuple{Int64},redFac;kargs...)

generates a Poisson disk sampling pattern for an Vector of size `shape` with a subsampling factor `redFac`.
The arguments are the same as in the 2d case with the exception that shape is of type `Tuple{Int64}`
"""
function sample_vdpoisson(shape::Tuple{Int64},redFac;kargs...)
  return sample_vdpoisson((shape[1],1),redFac;kargs...)
end

function sample_vdpoisson(shape::Tuple{Int64,Int64,Int64},redFac;kargs...)
  error("Not implemented")
end

"""
    sample_vdpoisson(shape::Tuple{Int64,Int64},redFac::Float64; rng=Random.default_rng(), seed=defaultSeed(rng),kargs...)

generates a variable density Poisson disk sampling pattern for an Array of size `shape` with a subsampling factor `redFac`.

# Arguments
* `shape::NTuple{2,Int64}` - size of the Array to be sampled
* `redFac::Float64`        - subsampling factor
* (`rng=Random.default_rng()`) - random number generator
* (`seed`)                 - seed for `rng`; with `nothing`, `rng` is used in its current state.
                             The default is 1234 for the default generator and `nothing` for any other
"""
function sample_vdpoisson(shape::Tuple{Int64,Int64},redFac::Float64; rng=Random.default_rng(), seed::Union{Nothing,Int64}=defaultSeed(rng),densityFac::Float64=1/pi,kargs...)
  M,N = shape
  samples = collect(1:M*N)
  chosenSamples = []
  forbiddenSamples =[]

  # set number of samples and mean distance
  numSamples = floor(Int64,M*N/redFac)
  rmean = densityFac*M*N/numSamples

  if !isnothing(seed)
    Random.seed!(rng, seed)
  end
  while length(chosenSamples)<numSamples
    # check if samples are available. Otherwise reduce rmean and recalculate samples
    if isempty(samples)
      rmean  = rmean/2.0
      forbiddenSamples = updateForbiddenSamples(M,N,chosenSamples,rmean)
      samples = setdiff(collect(1:M*N),forbiddenSamples)
      continue
    end

    # pick a sample from the available samples and add them to chosen samples
    chosen = samples[ rand(rng, 1:length(samples)) ]
    push!(chosenSamples,chosen)

    # calculate forbidden radius and points within
    fSphere = forbiddenSphere(M,N,chosen,rmean)
    forbiddenSamples = union(forbiddenSamples,fSphere)

    # update available samples
    samples = setdiff(samples,forbiddenSamples)
  end

  return chosenSamples
end

function rmin(kx::Float64, ky::Float64, ai::Float64, a0::Float64, rmean::Float64)
  m = rmean/(ai-3*a0)
  c = rmean*(ai-5*a0)/(2*(ai-3*a0))
  return m*sqrt(kx^2+ky^2)+c
end

function updateForbiddenSamples(M,N,samples,rmean)
  forbiddenSamples = []
  for i=1:length(samples)
    fSphere = forbiddenSphere(M,N,samples[i],rmean)
    forbiddenSamples = union(forbiddenSamples,fSphere)
  end
  return forbiddenSamples
end


function forbiddenSphere(M::Int64,N::Int64,chosen::Int64,rmean::Float64)
  x,y = Tuple(CartesianIndices((M,N))[chosen])
  kx = x-(M+1.)/2.0
  ky = y-(N+1.)/2.0
  rad_min = rmin(kx, ky, 0.5, 0.0, rmean)
  fsphere = []
  for xi = floor(Int64,x-M*rad_min):ceil(Int64,x+M*rad_min)
    # xi must be positive and smaller then M
    (xi<1) && continue
    (xi>M) && continue

    for yi =  floor(Int64,y-N*rad_min):ceil(Int64,y+N*rad_min)
      # yi must be positive and smaller then N
      (yi<1) && continue
      (yi>N) && continue

      if (xi-x)^2+(yi-y)^2 < rad_min^2
        idx = LinearIndices((M,N))[xi,yi]
        push!(fsphere,idx)
      end
    end
  end
  return fsphere
end

function generateMask(M,N,pat)
  mask = zeros(Int64,M,N)
  for i=1:length(pat)
    mask[pat[i]] = 1.0
  end
  return mask
end
