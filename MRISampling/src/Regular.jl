"""
    sample_regular(shape::Tuple{Int64,Int64}, redFac;kargs...)

generates a regular sampling pattern for an Array of size `shape` with a subsampling factor `redFac`.

# Arguments
* `shape::NTuple{N,Int64}` - size of the Array to be sampled
* `redFac::Float64`        - subsampling factor
"""
function sample_regular(shape::Tuple, redFac; kargs...)
  sample_regular(prod(shape), redFac; kargs...)
end

function sample_regular(maxInd::Int64, redFact; step=floor(Int,redFact), kargs...)
  outSet = range(1,step=step,length=floor(Int,maxInd/redFact))
  return outSet
end
