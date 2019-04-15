export WeightingOp

function WeightingOp(weights::Vector{T}, rep::Int=1) where T
  weights_cat = repeat(weights,rep)
  return opDiagonal(weights_cat)
end
