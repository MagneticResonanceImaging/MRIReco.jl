export SimplePatternParams

mutable struct SimplePatternParams
  nothing
end

function SimplePatternParams(;kargs...)
  SimplePatternParams(nothing)
end

function sample(shape::Tuple,redFac::Float64,patternParams::SimplePatternParams;kargs...)
  sample_simple(prod(shape),redFac)
end

function sample_simple(maxInd::Int64, redFact::Float64; kargs...)
  idxSet = collect(1:1:maxInd)
  outSet = zeros(Int64,floor(Int,maxInd/redFact))

  Random.seed!(1234)
  
  for i=1:floor(Int,maxInd/redFact)
    pick = rand(1:length(idxSet))
    outSet[i] = idxSet[pick]
    deleteat!(idxSet,pick)
  end

  return outSet
end
