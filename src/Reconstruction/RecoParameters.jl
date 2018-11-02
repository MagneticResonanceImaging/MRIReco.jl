export defaultRecoParams

function defaultRecoParams()
  params = Dict{Symbol,Any}()
  params[:reco] = "direct"
  params[:shape] = (32,32)
  params[:sparseTrafoName] = "Wavelet"
  params[:regularization] = "L1"
  params[:lambdL1] = 0.0
  params[:lambdL21] = 0.0
  params[:lambdTV] = 0.0
  params[:lambdLR] = 0.0
  params[:lambdLLR] = 0.0
  params[:normalizeReg] = false
  params[:solver] = "admm"
  params[:ρ] = 5.e-2
  params[:iterations] = 30

  return params
end
