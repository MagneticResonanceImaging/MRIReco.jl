

@testset "ESPIRiT" begin
  #image
  N=128
  img = shepp_logan(N)
  msk = zeros(N,N)
  msk[findall(x->x!=0,img)] .= 1

  # coil sensitivites
  smaps = birdcageSensitivity(N, 8, 1.5)
  snorm = sqrt.(sum(abs.(smaps).^2,dims=4))
  for i=1:8
    smaps[:,:,1,i] .= msk .* smaps[:,:,1,i] ./ snorm[:,:,1,1]
  end


  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N
  params[:senseMaps] = smaps

  acqData = simulation(img, params)
  acqData = MRIReco.sample_kspace(acqData, 2.0, "poisson", calsize=15)

  ksize=(6,6) # kernel size
  ncalib = 15 # number of calibration lines
  eigThresh_1 = 0.02  # threshold for picking singular vectors of calibration matrix
  eigThresh_2 = 0.95  # threshold for eigen vector decomposition in image space

  smaps2 = espirit(acqData,ksize,ncalib,eigThresh_1=eigThresh_1,eigThresh_2=eigThresh_2)

  # evaluate error only on the supprt of smaps
  for i=1:8
    smaps2[:,:,1,i] = msk .* smaps2[:,:,1,i]
  end
  # allow for a constant phase offset
  phs = mean(angle.(smaps))
  phs2 = mean(angle.(smaps2))
  smaps2 = exp(1im*(phs-phs2)) .* smaps2

  err = norm(vec(smaps2)-vec(smaps))/norm(vec(smaps))
  @test err < 3.e-2
end
