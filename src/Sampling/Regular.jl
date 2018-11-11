export RegularPatternParams

mutable struct RegularPatternParams
  nothing
end

function RegularPatternParams(;kargs...)
  SimplePatternParams(nothing)
end

function sample(shape::Tuple, redFac::Float64,
                patternParams::RegularPatternParams; kargs...)
  sample_regular(prod(shape), redFac)
end

function sample_regular(maxInd::Int64, redFact::Float64; kargs...)
  outSet = range(1,step=floor(Int,redFact),length=floor(Int,maxInd/redFact))
  return outSet
end
