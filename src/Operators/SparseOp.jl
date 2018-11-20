export SparseOp

#
# factory method for sparsifying operators
#
function SparseOp(name::AbstractString; kargs...)
  params = Dict(kargs)
  if name=="Wavelet"
    if get(params, :multiEcho, false)
      return MultiEchoWaveletOp(params[:shape], params[:numEchoes])
    else
      return WaveletOp(params[:shape])
    end
  elseif name=="nothing"
    return opEye(prod(params[:shape]))
  else
    error("SparseOp $name is not yet defined.")
  end

  return opEye(prod(params[:shape]))
end
