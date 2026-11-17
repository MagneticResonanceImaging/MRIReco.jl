export SparseOp

"""
    SparseOp(T::Type, name::AbstractString, shape::NTuple{N,Int64}; S = Vector{T}, kargs...) where N

generates the sparsifying transform (`<: AbstractLinearOperator`) given its name.

# Arguments
* `T::Type`                 - element type of the transform, for example `ComplexF64`
* `name::AbstractString`    - name of the sparsifying transform: `"Wavelet"` for `WaveletOp` with its default
                              wavelet, Daubechies 2, or `"nothing"` for the identity; other names throw an error
* `shape::NTuple{N,Int64}`  - size of the Array to be transformed
* (`S = Vector{T}`)         - storage vector type, for example `CuVector{ComplexF32}` on a GPU
* (`kargs`)                 - further keyword arguments, which are ignored
"""
function SparseOp(T::Type,name::AbstractString, shape::NTuple{N,Int64}; S = Vector{T}, kargs...) where N
  params = Dict(kargs)
  if name=="Wavelet"
    # if get(params, :multiEcho, false)
    #   return MultiEchoWaveletOp(params[:shape], params[:numEchoes])
    # else
    #   return WaveletOp(params[:shape])
    # end
    return WaveletOp(T; shape, S = S)
  elseif name=="nothing"
    return opEye(T,prod(shape), S = S)
  else
    error("SparseOp $name is not yet defined.")
  end

  return opEye(T, prod(shape), S = S)
end
