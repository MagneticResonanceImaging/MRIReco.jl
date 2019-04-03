export calculateIncoherence

function calculateIncoherence(acqData::AcquisitionData, recoParams::Dict, slice=1)

  N = prod(recoParams[:shape])

  @debug "Setup operator"
  F = EncodingOp2d(acqData, recoParams, slice=slice)
  @debug "Operator ready"

  maxValue = 0.0

  NRed = min(N, 500)

  pos = rand(1:N, NRed)

  @showprogress 1 "Computing Incoherence..." for i=1:length(pos)
    e = zeros(ComplexF64, N)
    e[pos[i]] = 1
    res = abs.( adjoint(F[1]) * (F[1] * e) )
    res ./= res[pos[i]]
    res[pos[i]] = 0
    maxValue = max(maxValue, maximum(res))
  end

  return maxValue
end
