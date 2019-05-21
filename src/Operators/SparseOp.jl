export SparseOp

"""
    SparseOp(name::AbstractString, shape::NTuple{N,Int64}; kargs...) where N

generates the sparsifying transform (`<: AbstractLinearOperator`) given its name.

# Arguments
* `name::AbstractString`    - name of the sparsifying transform
* `shape::NTuple{N,Int64}`  - size of the Array to be transformed
* (`kargs`)                 - additional keyword arguments
"""
function SparseOp(name::AbstractString, shape::NTuple{N,Int64}; kargs...) where N
  params = Dict(kargs)
  if name=="Wavelet"
    # if get(params, :multiEcho, false)
    #   return MultiEchoWaveletOp(params[:shape], params[:numEchoes])
    # else
    #   return WaveletOp(params[:shape])
    # end
    return WaveletOp(shape)
  elseif name=="nothing"
    return opEye(prod(shape))
  else
    error("SparseOp $name is not yet defined.")
  end

  return opEye(prod(shape))
end
