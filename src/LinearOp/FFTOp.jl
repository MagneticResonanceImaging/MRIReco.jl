export FFTOp

#
# Linear Operator to perform an FFT
#
function FFTOp(T::Type, shape::Tuple)
  plan = plan_fft(zeros(T, shape);flags=FFTW.MEASURE)
  iplan = plan_ifft(zeros(T, shape);flags=FFTW.MEASURE)

  return LinearOperator{T}(prod(shape), prod(shape), true, false
            , x->vec(fftshift(plan*fftshift(reshape(x,shape))))/sqrt(prod(shape))
            , nothing
            , y->vec(ifftshift(iplan*ifftshift(reshape(y,shape)))) * sqrt(prod(shape)) )
end
