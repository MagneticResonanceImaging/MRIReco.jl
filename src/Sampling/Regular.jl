function sample_regular(shape::Tuple, redFac::Float64; kargs...)
  sample_regular(prod(shape), redFac)
end

function sample_regular(maxInd::Int64, redFact::Float64; kargs...)
  outSet = range(1,step=floor(Int,redFact),length=floor(Int,maxInd/redFact))
  return outSet
end
