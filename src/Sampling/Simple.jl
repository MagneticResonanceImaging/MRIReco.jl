export SimplePatternParams

mutable struct SimplePatternParams
end

function SimplePatternParams(; kargs...)
  SimplePatternParams()
end

function sample(shape::Tuple,redFac::Float64, patternParams::SimplePatternParams;
                seed=1234,kargs...)
  sample_simple(prod(shape), redFac, seed)
end

function sample_simple(maxInd::Int64, redFact::Float64, seed)
  idxSet = collect(1:1:maxInd)
  outSet = zeros(Int64,floor(Int,maxInd/redFact))

  Random.seed!(seed)

  for i=1:floor(Int,maxInd/redFact)
    pick = rand(1:length(idxSet))
    outSet[i] = idxSet[pick]
    deleteat!(idxSet,pick)
  end

  return outSet
end
