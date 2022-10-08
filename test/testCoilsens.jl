

function testESPIRiT(N=128)

  #image
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

  # evaluate error only on the support of smaps
  for i=1:8
    smaps2[:,:,1,i] = msk .* smaps2[:,:,1,i]
  end

  # allow for a constant phase offset
  α = optimalScalingFactor(smaps, smaps2)
  # we don't want to allow scalings in the amplitude
  α /= abs(α)
  smaps2 .*= α

  err = norm(vec(smaps2)-vec(smaps))/norm(vec(smaps))
  @test err < 3.e-2

end


function testESPIRiT_newSize(imsize = 256)

  #normal image
  N=imsize÷2
  img = shepp_logan(N)
  msk = zeros(N,N)
  msk[findall(x->x!=0,img)] .= 1

  #larger image
  img2 = shepp_logan(imsize)
  msk2 = zeros(imsize,imsize)
  msk2[findall(x->x!=0,img2)] .= 1

  # coil sensitivites for normal image size
  smaps = birdcageSensitivity(N, 8, 1.5)
  snorm = sqrt.(sum(abs.(smaps).^2,dims=4))
  for i=1:8
    smaps[:,:,1,i] .= msk .* smaps[:,:,1,i] ./ snorm[:,:,1,1]
  end

  # coil sensitivites for larger image size
  smaps2 = birdcageSensitivity(imsize, 8, 1.5)
  snorm2 = sqrt.(sum(abs.(smaps2).^2,dims=4))
  for i=1:8
    smaps2[:,:,1,i] .= msk2 .* smaps2[:,:,1,i] ./ snorm2[:,:,1,1]
  end
  
  # simulation for normal image size
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

  emaps = espirit(acqData,ksize,ncalib,(imsize,imsize),eigThresh_1=eigThresh_1,eigThresh_2=eigThresh_2)

  # simulation for larger image size
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, imsize)
  params[:numSamplingPerProfile] = imsize
  params[:senseMaps] = smaps2

  acqData2 = simulation(img2, params)
  acqData2 = MRIReco.sample_kspace(acqData2, 2.0, "poisson", calsize=15)
  emaps2 = espirit(acqData2,ksize,ncalib,(imsize,imsize),eigThresh_1=eigThresh_1,eigThresh_2=eigThresh_2)

  # evaluate error only on the support of smaps
  for i=1:8
    emaps2[:,:,1,i] = msk2 .* emaps2[:,:,1,i]
    emaps[:,:,1,i] = msk2 .* emaps[:,:,1,i]
  end

  # allow for a constant phase offset
  α = optimalScalingFactor(emaps, emaps2)
  # we don't want to allow scalings in the amplitude
  α /= abs(α)
  emaps2 .*= α

  err = norm(vec(emaps2)-vec(emaps))/norm(vec(emaps))
  @test err < 3.e-2


  ## estimateCoilSensitivities
  acqData = simulation(img, params)
  params[:reco] = "direct"
  params[:reconSize] = (N,N)
  Ireco = reconstruction(acqData,params)

  smaps3 = estimateCoilSensitivities(Ireco)

  msk = abs.(smaps3[:,:,1,1,1,1]).>0

  err = nrmsd(smaps .* msk, msk .* smaps3[:,:,:,1,:,1])
  @test err < 1e-15
end

@testset "ESPIRiT" begin

  # test default ESPIRiT using acqData.encodingSize as output map size
  testESPIRiT()

  # test making maps of new size from low-res acq data
  # even matrix size
  testESPIRiT_newSize(128)

  # odd matrix size
  testESPIRiT_newSize(127)

end
