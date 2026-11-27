"""
    sample_vdpoisson_bart(shape::Tuple{Int64},redFac;kargs...)

generates a variable density Poisson disk sampling pattern, as the `poisson` command of BART, for a Vector of size
`shape` with a subsampling factor `redFac`.
The arguments are the same as in the 2d case with the exception that shape is of type `Tuple{Int64}`
"""
function sample_vdpoisson_bart(shape::Tuple{Int64},redFac;kargs...)
  return sample_vdpoisson_bart((shape[1],1),redFac;kargs...)
end

function sample_vdpoisson_bart(shape::Tuple{Int64,Int64,Int64},redFac;kargs...)
  error("Not implemented")
end

"""
    sample_vdpoisson_bart(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0,vardens::Real=20.0,
                          elliptical::Bool=false,rng=Random.default_rng(),seed=nothing,kargs...)

generates a variable density Poisson disk sampling pattern for an Array of size `shape` with a subsampling factor
`redFac`, with the density law of the `poisson` command of the Berkeley Advanced Reconstruction Toolbox (BART).

The points are placed with Bridson's algorithm in a continuous generation space: each new point is proposed in the
annulus between one and two times the minimum distance of an existing point, and accepted if no point lies within its
own minimum distance

    d(p) = δ (1 + vardens |p - c|²),   δ = 1 / (1.275 max(M, N)),

where `c` is the center of the generation space. The generation space is then stretched by an acceleration factor `a`
along both dimensions (the options `-y a -z a` of BART), and the points are rounded down to the grid. In the final
pattern, the minimum distance is therefore `0.784 (a + vardens ρ² / (4a))` grid points, with `ρ` the distance to the
center of k-space normalized to 1 at the edge. `a` is adjusted so that the pattern contains `floor(M*N/redFac)`
points, the calibration area included; if the adjusted pattern has a few points too many, randomly chosen points
outside the calibration area are removed. The center of k-space (index `(M÷2+1, N÷2+1)`) is always sampled, also
with `calsize = 0`.

Because `vardens` acts in the generation space, which `a` stretches, the number of points is not monotonic in `a`,
and each `vardens` has a smallest reachable reduction factor (about 8 for `vardens = 20` and `elliptical = true`).
A smaller `redFac` throws an `ArgumentError`; use a smaller `vardens` in this case.

# Arguments
* `shape::NTuple{2,Int64}`    - size of the Array to be sampled
* `redFac::Float64`           - subsampling factor
* (`calsize::Int64=0`)        - size of the fully sampled calibration area
* (`vardens::Real=20.0`)      - strength of the variable density (option `-V` of BART; `-v` sets 20; 0 gives a
                                uniform Poisson disk pattern)
* (`elliptical::Bool=false`)  - only sample inside the ellipse inscribed in the grid (option `-e` of BART)
* (`rng=Random.default_rng()`) - random number generator
* (`seed=nothing`)            - seed for `rng`; with `nothing`, `rng` is used in its current state

# References
* BART, `src/poisson.c` and `src/misc/pd.c`: https://github.com/mrirecon/bart
* R. Bridson, Fast Poisson disk sampling in arbitrary dimensions, SIGGRAPH sketches, 2007.
"""
function sample_vdpoisson_bart(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0,vardens::Real=20.0,
                               elliptical::Bool=false,rng=Random.default_rng(),
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
  pattern(a) = bartpoisson_pattern(M,N,a,Float64(vardens),elliptical,calib,MersenneTwister(s))

  # densest pattern: coarse search over a
  as = exp.(range(log(0.5), log(3.0*max(1.0,sqrt(vardens/20))), length=7))
  pats = [pattern(a) for a in as]
  i = argmax(count.(pats))
  if count(pats[i]) < NumberSamples
    throw(ArgumentError("the reduction factor $redFac cannot be reached with vardens = $vardens; the smallest " *
                        "reachable reduction factor is about $(round(M*N/count(pats[i]), digits=2))"))
  end

  # the number of points decreases with a beyond the densest pattern: bisection
  lo, hi = as[i], 2*as[i]
  P = pattern(hi)
  while count(P) >= NumberSamples
    lo, hi = hi, 2*hi
    P = pattern(hi)
  end
  best = pattern(lo)
  for _ = 1:25
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
Bridson's algorithm with the density law of BART, stretched by `a` and rounded down to an `M` by `N` grid.
`calib` marks the calibration area, which is always sampled.
"""
function bartpoisson_pattern(M::Int64,N::Int64,a::Float64,vardens::Float64,elliptical::Bool,calib::BitMatrix,
                             rng::AbstractRNG; k::Int64=30)
  K = max(M,N)
  δ = (1/1.275)/K
  sy, sz = a*K/M, a*K/N
  # only the part of the generation space that lands in the grid is generated
  lo = (0.5-0.5/sy, 0.5-0.5/sz)
  hi = (0.5+0.5/sy, 0.5+0.5/sz)
  scale(p) = 1 + vardens*((p[1]-0.5)^2 + (p[2]-0.5)^2)

  # background grid with at most one point per cell
  cell = δ/sqrt(2)
  G = (ceil(Int64,(hi[1]-lo[1])/cell)+1, ceil(Int64,(hi[2]-lo[2])/cell)+1)
  grid = zeros(Int64,G)
  gpos(p) = (floor(Int64,(p[1]-lo[1])/cell)+1, floor(Int64,(p[2]-lo[2])/cell)+1)

  pts = [(0.5,0.5)]
  grid[gpos(pts[1])...] = 1
  active = [1]
  while !isempty(active)
    ia = rand(rng, 1:length(active))
    s = pts[active[ia]]
    dd = δ*scale(s)
    found = false
    for _ = 1:k
      p = s
      while true # candidate in the annulus [dd, 2dd]
        p = (s[1]+(rand(rng)-0.5)*4*dd, s[2]+(rand(rng)-0.5)*4*dd)
        (lo[1] <= p[1] <= hi[1] && lo[2] <= p[2] <= hi[2]) || continue
        dd <= hypot(p[1]-s[1], p[2]-s[2]) <= 2*dd && break
      end
      r = δ*scale(p)
      c = gpos(p)
      n = ceil(Int64,r/cell)
      accept = true
      for j = max(1,c[2]-n):min(G[2],c[2]+n), i = max(1,c[1]-n):min(G[1],c[1]+n)
        q = grid[i,j]
        if q != 0 && hypot(pts[q][1]-p[1], pts[q][2]-p[2]) <= r
          accept = false
          break
        end
      end
      accept || continue
      push!(pts,p)
      grid[c...] = length(pts)
      push!(active,length(pts))
      found = true
      break
    end
    found || deleteat!(active,ia)
  end

  mask = copy(calib)
  for p in pts
    q = ((p[1]-0.5)*sy+0.5, (p[2]-0.5)*sz+0.5)
    d = elliptical ? hypot(q[1]-0.5, q[2]-0.5) : max(abs(q[1]-0.5), abs(q[2]-0.5))
    d <= 0.5 || continue
    y, z = floor(Int64,q[1]*M)+1, floor(Int64,q[2]*N)+1
    if 1 <= y <= M && 1 <= z <= N
      mask[y,z] = true
    end
  end
  return mask
end

"""
removes randomly chosen samples outside the calibration area `calib` from `mask`
until it contains `NumberSamples` samples
"""
function removeexcesssamples!(mask::BitMatrix,calib::BitMatrix,NumberSamples::Int64,rng::AbstractRNG)
  excess = count(mask) - NumberSamples
  if excess > 0
    candidates = findall(mask .& .!calib)
    mask[candidates[randperm(rng,length(candidates))[1:excess]]] .= false
  end
  return mask
end
