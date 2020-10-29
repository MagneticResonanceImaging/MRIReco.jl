using StatsBase

"""
  use the Loyd-Max algorithm to find a quantizer with `M` quantization levels
  for the signal `sig`.
"""
function loydMax(sig::Vector{T}, M::Int64; iterations=100) where T<:Real
  # sort signal to simplify the algorithm
  sig_sort = sort(sig)

  # start with uniform levels
  x = collect(range(minimum(sig), stop=maximum(sig),length=M))
  for i=1:iterations
    # compute decision thresholds
    t = 0.5*(x[1:end-1].+x[2:end])
    # @info "t = $(t)"
    # update representative levels
    num = 0 
    denom = 0
    idx = 1
    for j=1:length(sig)
      if idx<M && sig_sort[j]<t[idx]
        num += sig_sort[j]
        denom += 1
      elseif idx<M
        x[idx] = num/denom
        num = 0
        denom = 0
        idx += 1
      else
        num += sig_sort[j]
        denom += 1
      end
      x[M] = num/denom
    end
  end

  return x
end