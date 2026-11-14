function sample_random(shape::Tuple{Int64},redFac;calsize::Int64=0,kargs...)
  return sample_random((shape[1],1),redFac;calsize=calsize,kargs...)
end

function sample_random(shape::Tuple{Int64,Int64,Int64},redFac;calsize::Int64=0,kargs...)
  error("Not implemented")
end

"""
    sample_random(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0,rng=Random.default_rng(),seed=nothing,kargs...)

generates a random sampling pattern for an Array of size `shape` with a subsampling factor `redFac`.

# Arguments
* `shape::NTuple{N,Int64}` - size of the Array to be sampled
* `redFac::Float64`        - subsampling factor
* (`calsize::Int64=0`)     - size of the fully sampled calibration area
* (`rng=Random.default_rng()`) - random number generator
* (`seed=nothing`)         - seed for `rng`; with `nothing`, `rng` is used in its current state
"""
function sample_random(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0,rng=Random.default_rng(),seed::Union{Nothing,Int64}=nothing,kargs...)
  M,N = shape
  A = zeros(Int64,M,N)
  ChosenSamples = 0
  NumberSamples = floor(Int64,M*N/redFac)
  A, ChosenSamples = setcalibrationarea(A,M,N,redFac,calsize)

  chosen = (LinearIndices(A))[findall(x->x!=0, A)]
  B = abs.(A .- 1)
  selection = (LinearIndices(B))[findall(x->x!=0, B)]
  if !isnothing(seed)
    Random.seed!(rng, seed)
  end
  while ChosenSamples < NumberSamples
    randnum = rand(rng, 1:length(selection))
    push!(chosen,selection[randnum])
    deleteat!(selection,randnum)
    ChosenSamples += 1
  end
  A = reshape(A,M*N,1)
  A[chosen] .= 1

  A = reshape(A,M,N)

  return chosen
end
