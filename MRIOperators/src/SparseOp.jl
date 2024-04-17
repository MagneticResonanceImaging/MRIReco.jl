export SparseOp

"""
    SparseOp(type::T,name::AbstractString, shape::NTuple{N,Int64}; kargs...) where{N,T}

generates the sparsifying transform (`<: AbstractLinearOperator`) given its name.

# Arguments
* `name::AbstractString`    - name of the sparsifying transform
* `shape::NTuple{D,Int64}`  - size of the Array to be transformed
* (`kargs`)                 - additional keyword arguments
"""
function SparseOp(T::Type,name::AbstractString, shape::NTuple{N,Int64}; S = Vector{T}, kargs...) where N
  params = Dict(kargs)
  if name=="Wavelet"
    # if get(params, :multiEcho, false)
    #   return MultiEchoWaveletOp(params[:shape], params[:numEchoes])
    # else
    #   return WaveletOp(params[:shape])
    # end
    return WaveletOp(T; shape) # does not support S
  elseif name=="nothing"
    return opEye(T,prod(shape), S = S)
  else
    error("SparseOp $name is not yet defined.")
  end

  return opEye(T, prod(shape), S = S)
end
