function setcalibrationarea(A::Array{Int64,2},M::Int64,N::Int64,redFac::Float64,calsize::Int64)
  cal_x, cal_y = min(calsize,M), min(calsize,N)
  ChosenSamples = 0
  NumberSamples = floor(Int64,M*N/redFac)
  if calsize > 0
    if cal_x*cal_y > NumberSamples
      @error "Number of calibration samples exceeds number of required total samples"
    end
    ChosenSamples = cal_x*cal_y
    A[floor(Int64,end/2)-floor(Int64,cal_x/2)+1:floor(Int64,end/2)-floor(Int64,cal_x/2)+cal_x,floor(Int64,end/2)-floor(Int64,cal_y/2)+1:floor(Int64,end/2)-floor(Int64,cal_y/2)+cal_y].=1
  end

  return A, ChosenSamples

end
