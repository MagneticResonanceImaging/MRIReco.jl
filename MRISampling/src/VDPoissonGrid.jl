"""
    sample_vdpoisson_grid(shape::Tuple{Int64},redFac;kargs...)

generates a variable density Poisson disk sampling pattern on the grid for a Vector of size `shape` with a
subsampling factor `redFac`.
The arguments are the same as in the 2d case with the exception that shape is of type `Tuple{Int64}`
"""
function sample_vdpoisson_grid(shape::Tuple{Int64},redFac;kargs...)
  return sample_vdpoisson_grid((shape[1],1),redFac;kargs...)
end

function sample_vdpoisson_grid(shape::Tuple{Int64,Int64,Int64},redFac;kargs...)
  error("Not implemented")
end

"""
    sample_vdpoisson_grid(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0,vardens::Real=2.0,jitter::Real=1.0,
                          elliptical::Bool=false,rng=Random.default_rng(),seed=nothing,kargs...)

generates a variable density Poisson disk sampling pattern for an Array of size `shape` with a subsampling factor
`redFac`, by dart throwing on the grid.

All points of the grid are visited once in a random order. A point is kept if no kept point lies closer than

    r(ρ) = r0 (1 + vardens ρ),

where `ρ` is the distance to the center of k-space normalized to 1 at the edge. The distances are computed with
positions that are shifted by a random offset of up to ± `jitter`/2 grid points: on the integer grid, only a few
distances exist, and without this offset a variable density gives regular rings and checkerboard patterns. `r0` is
adjusted so that the pattern contains `floor(M*N/redFac)` points, the calibration area included; if the adjusted
pattern has a few points too many, randomly chosen points outside the calibration area are removed. The center of
k-space (index `(M÷2+1, N÷2+1)`) is always sampled, also with `calsize = 0`. The points of the calibration area and
the center do not exclude their neighbours.

A large `vardens` does not allow small reduction factors, because the edge cannot be denser than the grid; a smaller
`redFac` throws an `ArgumentError`.

# Arguments
* `shape::NTuple{2,Int64}`    - size of the Array to be sampled
* `redFac::Float64`           - subsampling factor
* (`calsize::Int64=0`)        - size of the fully sampled calibration area
* (`vardens::Real=2.0`)       - strength of the variable density, the slope of the linear law (0 gives a uniform
                                Poisson disk pattern)
* (`jitter::Real=1.0`)        - random offset of the positions for the distances, in grid points
* (`elliptical::Bool=false`)  - only sample inside the ellipse inscribed in the grid
* (`rng=Random.default_rng()`) - random number generator
* (`seed=nothing`)            - seed for `rng`; with `nothing`, `rng` is used in its current state
"""
function sample_vdpoisson_grid(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0,vardens::Real=2.0,
                               jitter::Real=1.0,elliptical::Bool=false,rng=Random.default_rng(),
                               seed::Union{Nothing,Int64}=nothing,kargs...)
  M,N = shape
  NumberSamples = floor(Int64,M*N/redFac)
  A, _ = setcalibrationarea(zeros(Int64,M,N),M,N,redFac,calsize)
  calib = A .!= 0
  calib[M÷2+1,N÷2+1] = true # the center of k-space is always sampled, also without calibration area

  if !isnothing(seed)
    Random.seed!(rng, seed)
  end
  # all the trials of the search use the same random numbers; rng advances by one draw
  s = rand(rng, UInt64)
  pattern(r0) = gridpoisson_pattern(M,N,r0,Float64(vardens),Float64(jitter),elliptical,calib,MersenneTwister(s))

  # the number of points decreases with r0: bisection
  lo = 0.3
  best = pattern(lo)
  if count(best) < NumberSamples
    throw(ArgumentError("the reduction factor $redFac cannot be reached with vardens = $vardens; the smallest " *
                        "reachable reduction factor is about $(round(M*N/count(best), digits=2))"))
  end
  hi = 2*sqrt(M*N/NumberSamples)
  while count(pattern(hi)) >= NumberSamples
    lo, hi = hi, 2*hi
  end
  for _ = 1:40
    mid = (lo+hi)/2
    P = pattern(mid)
    if count(P) >= NumberSamples
      lo, best = mid, P
    else
      hi = mid
    end
    count(best) - NumberSamples <= max(1, NumberSamples ÷ 1000) && break
  end

  removeexcesssamples!(best,calib,NumberSamples,MersenneTwister(s+1))
  return (LinearIndices(best))[findall(best)]
end

"""
dart throwing on an `M` by `N` grid with the minimum distance `r0 (1 + vardens ρ)`; `calib` marks the calibration area,
which is always sampled.
"""
function gridpoisson_pattern(M::Int64,N::Int64,r0::Float64,vardens::Float64,jitter::Float64,elliptical::Bool,
                             calib::BitMatrix,rng::AbstractRNG)
  jy = jitter .* (rand(rng,M,N) .- 0.5)
  jz = jitter .* (rand(rng,M,N) .- 0.5)
  cy, cz = M÷2+1, N÷2+1
  mask = copy(calib)
  kept = falses(M,N)
  for idx in randperm(rng,M*N)
    y, z = Tuple(CartesianIndices((M,N))[idx])
    mask[y,z] && continue
    ρ = sqrt(((y-cy)/(M/2))^2 + ((z-cz)/(N/2))^2)
    (elliptical && ρ > 1) && continue
    r = r0*(1 + vardens*ρ)
    d = ceil(Int64,r+jitter)
    accept = true
    for dz = -d:d, dy = -d:d
      yy, zz = y+dy, z+dz
      (1 <= yy <= M && 1 <= zz <= N) || continue
      if kept[yy,zz] && (dy+jy[yy,zz]-jy[y,z])^2 + (dz+jz[yy,zz]-jz[y,z])^2 < r^2
        accept = false
        break
      end
    end
    accept || continue
    kept[y,z] = true
    mask[y,z] = true
  end
  return mask
end
