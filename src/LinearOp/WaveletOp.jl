export WaveletOp, MultiEchoWaveletOp

#
# Linear Operator to perform Wavelet transformation
#
# function WaveletOp(shape, wt=wavelet(WT.db2))
#   return LinearOperator{ComplexF64}(prod(shape), prod(shape), false, false
#             , x->vec( dwt(reshape(x, shape), wt) )
#             , Nullable{Function}()
#             , y->vec( idwt( reshape(y, shape), wt) ) )
# end

function WaveletOp(shape, wt=wavelet(WT.db2))
  return LinearOperator{ComplexF64,Function,Nothing,Function}(maximum(shape)^2, prod(shape), false, false
            , x->waveletProd(x,shape,wt)
            , nothing
            , y->waveletCTProd(y,shape,wt) )
end

function MultiEchoWaveletOp(shape, numEchoes, wt=wavelet(WT.db2))
  W = WaveletOp(shape, wt)
  return diagOp([W for i=1:numEchoes]...)
end

function waveletProd(x::Vector{T},shape, wt) where T
  if shape[1] != shape[2]
    xSquare = zeros(T, maximum(shape), maximum(shape))
    xSquare[1:shape[1],1:shape[2]] = reshape(x,shape)
  else
    xSquare = reshape(x,shape)
  end
  return vec( dwt(xSquare, wt) )
end

function waveletCTProd(y::Vector{T},shape, wt) where T
  squareSize = (maximum(shape), maximum(shape))
  x = idwt( reshape(y, squareSize), wt)
  return vec( x[1:shape[1],1:shape[2]] )
end
