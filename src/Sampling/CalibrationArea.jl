function setcalibrationarea(A::Array{Int64,2},M::Int64,N::Int64,redFac::Float64,calsize::Int64)
  if calsize > M || calsize > N
    error("Calibration size larger than dimension of sampling matrix")
  end
  ChosenSamples = 0
  NumberSamples = floor(Int64,M*N/redFac)
  if calsize > 0
    if calsize*calsize > NumberSamples
      error("Number of calibration samples exceeds number of required total samples")
    end
    ChosenSamples = calsize*calsize
    A[floor(Int64,end/2)-floor(Int64,calsize/2)+1:floor(Int64,end/2)-floor(Int64,calsize/2)+calsize,floor(Int64,end/2)-floor(Int64,calsize/2)+1:floor(Int64,end/2)-floor(Int64,calsize/2)+calsize].=1
  end

  return A, ChosenSamples

end
