export SparseOp

#
# factory method for sparsifying operators
#
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
