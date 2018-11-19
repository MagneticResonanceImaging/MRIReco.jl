export LinesPatternParams

mutable struct LinesPatternParams
  nothing
end

function LinesPatternParams(;kargs...)
  LinesPatternParams(nothing)
end

function sample(shape::Tuple{Int64},redFac,patternParams::LinesPatternParams;kargs...)
  error("Not implemented")
end

function sample(shape::Tuple{Int64,Int64},redFac,patternParams::LinesPatternParams;kargs...)
  sample_lines(shape[1],shape[2],redFac;kargs...)
end

function sample(shape::Tuple{Int64,Int64,Int64},redFac,patternParams::LinesPatternParams;kargs...)
  sample_lines_3D(shape[1],shape[2],shape[3],redFac;kargs...)
end

function sample_lines(M::Int64,N::Int64,redFac::Float64;sampleFunc="simple",kargs...)

  A = zeros(Int64,M,N)

  #yInd = sample((M,),redFac,SimplePatternParams())

  # OLD version with exchanged order of phase and frequency encoding
  # if sampleFunc=="simple"
  #   yInd = sample((M,),redFac,SimplePatternParams())
  # elseif sampleFunc =="vardens"
  #   yInd = sample((M,),redFac,VardensPatternParams(;kargs...))
  # elseif sampleFunc =="poisson"
  #   yInd = sample((M,),redFac,PoissonDiskPatternParams(;kargs...))
  # else
  #   error("No valid sampling pattern specified")
  # end
  #
  # A[yInd,1:N] = 1;

  if sampleFunc=="simple"
    yInd = sample((N,),redFac,SimplePatternParams(;kargs...); kargs...)
  elseif sampleFunc=="regular"
    yInd = sample((N,),redFac,RegularPatternParams())
  elseif sampleFunc =="vardens"
    yInd = sample((N,),redFac,VardensPatternParams(;kargs...))
  elseif sampleFunc =="poisson"
    yInd = sample((N,),redFac,PoissonDiskPatternParams(;kargs...))
  else
    error("No valid sampling pattern specified")
  end

  A[1:M,yInd] .= 1;

  return (LinearIndices(A))[findall(x->x!=0, A)] # find(A)

end

function sample_lines_3D(M::Int64,N::Int64,Z::Int64,redFac::Float64;xyFunc="simple",kargs...)
  A2D = zeros(Int64,M,Z)
  A3D = zeros(Int64,M,N,Z)

  if xyFunc=="simple"
    idx = sample(size(A2D),redFac,SimplePatternParams())
  elseif xyFunc=="regular"
    idx = sample(size(A2D),redFac,RegularPatternParams())
  elseif xyFunc =="vardens"
    idx = sample(size(A2D),redFac,VardensPatternParams(;kargs...))
  elseif xyFunc =="poisson"
    idx = sample(size(A2D),redFac,PoissonDiskPatternParams(;kargs...))
  else
    error("No valid 2D sampling pattern specified")
  end


  A2D[idx] = 1;

  for i=1:N
    A3D[:,i,:] = A2D[:]
  end

  return find(A3D)

end
