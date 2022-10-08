export WeightingOp

"""
    WeightingOp(weights::Vector{T}, rep::Int=1) where T

generates a `LinearOperator` which multiplies an input vector `weights`

# Arguments
* `weights::Vector{T}` - weights vector
* `rep::Int=1`         - number of sub-arrays that need to be multiplied with `weights`
"""
function WeightingOp(weights::Vector{T}, rep::Int=1) where T
  weights_cat = repeat(weights,rep)
  return opDiagonal(weights_cat)
end
