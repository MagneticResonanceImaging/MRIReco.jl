function sample_lines(shape::Tuple{Int64},redFac;sampleFunc="random",kargs...)
  error("Not implemented")
end

function sample_lines(shape::Tuple{Int64,Int64},redFac::Float64;sampleFunc="random",kargs...)
  M,N=shape
  A = zeros(Int64,shape)
  yInd = sample((N,),redFac,sampleFunc;kargs...)
  A[1:M,yInd] .= 1;

  return (LinearIndices(A))[findall(x->x!=0, A)]
end

function sample_lines(shape::Tuple{Int64,Int64,Int64},redFac::Float64;sampleFunc="random",kargs...)
  M,N,Z=shape
  A2D = zeros(Int64,N,Z)
  A3D = zeros(Int64,shape)

  idx = sample(size(A2D),redFac,sampleFunc;kargs...)

  A2D[idx] .= 1;

  for i=1:M
    A3D[i,:,:] = A2D[:]
  end

  return (LinearIndices(A3D))[findall(x->x!=0, A3D)]
end
