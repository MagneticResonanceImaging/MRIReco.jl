export VardensPatternParams

mutable struct VardensPatternParams
  calsize::Int64
end

function VardensPatternParams(;calsize::Int64=0,kargs...)
  VardensPatternParams(calsize)
end

function sample(shape::Tuple{Int64},redFac,patternParams::VardensPatternParams;kargs...)
  error("Not implemented")
end

function sample(shape::Tuple{Int64,Int64},redFac,patternParams::VardensPatternParams;kargs...)
  sample_vardens(shape[1],shape[2],redFac,calsize=patternParams.calsize)
end

function sample(shape::Tuple{Int64,Int64,Int64},redFac,patternParams::VardensPatternParams;kargs...)
  error("Not implemented")
end

function sample_vardens(M::Int64,N::Int64,redFac::Float64;calsize::Int64=0,kargs...)

  A = zeros(Int64,M,N)
  ChosenSamples = 0
  NumberSamples = floor(Int64,M*N/redFac)
  A, ChosenSamples = setcalibrationarea(A,M,N,redFac,calsize)

  chosen = (LinearIndices(A))[findall(x->x!=0, A)]
  B = abs.(A .- 1)
  selection = (LinearIndices(B))[findall(x->x!=0, B)]
  #srand(1234)
  while ChosenSamples < NumberSamples
    randnum = rand(1:length(selection))
    push!(chosen,selection[randnum])
    deleteat!(selection,randnum)
    ChosenSamples += 1
  end
  A = reshape(A,M*N,1)
  A[chosen] .= 1

  A = reshape(A,M,N)

  return chosen
end
