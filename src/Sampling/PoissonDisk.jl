export PoissonDiskPatternParams

mutable struct PoissonDiskPatternParams
  calsize::Int64
  seed::Int64
end

function PoissonDiskPatternParams(;calsize::Int64=0, seed::Int64 = 1234, kargs...)
  PoissonDiskPatternParams(calsize, seed)
end

function sample(shape::Tuple{Int64},redFac,patternParams::PoissonDiskPatternParams;kargs...)
  sample_poissondisk(shape[1],one(Int64),redFac,calsize=patternParams.calsize,seed=patternParams.seed;kargs...)
end

function sample(shape::Tuple{Int64,Int64},redFac,patternParams::PoissonDiskPatternParams;kargs...)
  sample_poissondisk(shape[1],shape[2],redFac,calsize=patternParams.calsize,seed=patternParams.seed;kargs...)
end

function sample(shape::Tuple{Int64,Int64,Int64},redFac,patternParams::PoissonDiskPatternParams;kargs...)
  error("Not implemented")
end

function sample_poissondisk(M::Int64,N::Int64,redFac::Float64;calsize::Int64=0, seed::Int64=1234,kargs...)
  A = zeros(Int64,M,N)
  ChosenSamples = 0
  NumberSamples = floor(Int64,M*N/redFac)
  A, ChosenSamples = setcalibrationarea(A,M,N,redFac,calsize)
  chosen = (LinearIndices(A))[findall(x->x!=0, A)]
  # selection = find(abs.(A .- 1))
  selection = (LinearIndices(abs.(A .- 1)))[findall(x->x!=0, abs.(A .- 1))]

  forbiddenradius = floor(Int64,sqrt((M*N-calsize^2)/(NumberSamples-ChosenSamples)/pi))
  if forbiddenradius <= 1
    forbiddenradius = 2
  end
 #println(forbiddenradius)

  Random.seed!(seed)

  while ChosenSamples < NumberSamples
    if length(selection)<1
      A = reshape(A,M*N,1)
      A[chosen] .= 1
      A = reshape(A,M,N)
      selection = (LinearIndices(abs.(A .- 1)))[findall(x->x!=0, abs.(A .- 1))]
    end
    randnum = rand(1:length(selection))
    push!(chosen,selection[randnum])
    ChosenSamples += 1

    selection = deleteinsideforbiddenradius(chosen[end],selection,forbiddenradius,M,N)
  end
  A = reshape(A,M*N,1)
  A[chosen] .= 1

  A = reshape(A,M,N)

  return chosen
end

function deleteinsideforbiddenradius(chosen::Int64,selection::Array{Int64},rad::Int64,M::Int64,N::Int64)
  x,y=Tuple(CartesianIndices((M,N))[chosen])

  for xi = x-rad:x+rad
      for yi = y-rad:y+rad
        if (xi > 0 && xi<= M && yi > 0 && yi<= N)
        index = (LinearIndices((M,N)))[xi,yi]
        deleteindex = 0
        for (indexi,valuei) in enumerate(selection)
           if valuei == index
              deleteindex = indexi
           end
        end
        if deleteindex > 0
          deleteat!(selection,deleteindex)
        end
      end
      end
  end
  return selection
end
