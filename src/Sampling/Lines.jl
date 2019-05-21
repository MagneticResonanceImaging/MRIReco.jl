function sample_lines(shape::Tuple{Int64},redFac;sampleFunc="random",kargs...)
  error("Not implemented")
end

"""
    sample_lines(shape::Tuple{Int64,Int64},redFac::Float64;sampleFunc="random",kargs...)

generates a pattern to sample complete lines of an Array of size `shape` with a subsampling factor `redFac`.

# Arguments
* `shape::NTuple{N,Int64}` - size of the Array to be sampled
* `redFac::Float64`        - subsampling factor
* `sampleFunc="random"`        - name of the sampling function
                            ("random, "regular", "lines", "poisson" or "vdPoisson")
* `kargs...`               - addional keyword arguments
"""
function sample_lines(shape::Tuple{Int64,Int64},redFac::Float64;sampleFunc="random",kargs...)
  M,N=shape
  A = zeros(Int64,shape)
  yInd = sample((N,),redFac,sampleFunc;kargs...)
  A[1:M,yInd] .= 1;

  return (LinearIndices(A))[findall(x->x!=0, A)]
end

"""
    sample_lines(shape::Tuple{Int64,Int64,Int64},redFac::Float64;sampleFunc="random",kargs...)

generates a pattern to sample complete lines of an Array of size `shape` with a subsampling factor `redFac`.
The arguments are the same as in the 2d case with the exception that shape is of type `NTuple{3,Int64}`
"""
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
