export defaultRecoParams

function defaultRecoParams()
  params = Dict{Symbol,Any}()
  params[:reco] = "direct"
  params[:reconSize] = (32,32)
  params[:sparseTrafoName] = "Wavelet"
  params[:regularization] = "L1"
  params[:λ] = 0.0
  params[:normalizeReg] = false
  params[:solver] = "admm"
  params[:ρ] = 5.e-2
  params[:iterations] = 30

  return params
end
