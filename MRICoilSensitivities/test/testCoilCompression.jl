function test2DCC(N=64)
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

  #reconstruction
  params = Dict{Symbol,Any}()
  params[:reco] = "direct"
  params[:reconSize] = acqData.encodingSize

  x_ori = reconstruction(acqData, params)
  x_ori = mergeChannels(x_ori)

  acqDataCC,ccMat = softwareCoilCompression(acqData, 6)
  x_cc = reconstruction(acqDataCC, params)
  x_cc = mergeChannels(x_cc)

  @test (norm(vec(x_ori)-vec(x_cc))/norm(vec(x_ori))) < 2e-4
end

function test3DCC(N=64, NSl=64)
  ## 3D espirit
  type = ComplexF32
  sh = shepp_logan(N)
  x = repeat(sh,1,1,NSl)
  msk = zeros(type,N,N,NSl)
  msk[findall(x->x!=0,x)] .= 1

  # coil sensitivites
  smaps = birdcageSensitivity(N, 8, 1.5)
  smaps = repeat(smaps,1,1,NSl,1)
  snorm = sqrt.(sum(abs.(smaps).^2,dims=4))
  for i=1:8
      smaps[:,:,:,i] .= msk .* smaps[:,:,:,i] ./ snorm[:,:,:,1]
  end

  # convert to type
  x = convert.(type,x)
  smaps = convert.(type,smaps)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian3D"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N
  params[:numSlices] = NSl
  params[:senseMaps] = smaps

  acqData = simulation( x, params )

  #reconstruction
  params = Dict{Symbol,Any}()
  params[:reco] = "direct"
  params[:reconSize] = acqData.encodingSize

  x_ori = reconstruction(acqData, params)
  x_ori = mergeChannels(x_ori)

  acqDataCC,ccMat = softwareCoilCompression(acqData, 6)
  x_cc = reconstruction(acqDataCC, params)
  x_cc = mergeChannels(x_cc)

  @test (norm(vec(x_ori)-vec(x_cc))/norm(vec(x_ori))) < 2e-4

  # extract kspace then perfom Coil Compression
   kdata = kDataCart(acqData)

   kdataCC,ccMat2 = softwareCoilCompression(kdata, 6)
   acqDataCC2 = AcquisitionData(kdataCC)
   x_cc2 = reconstruction(acqDataCC2, params)
   x_cc2 = mergeChannels(x_cc2)
   @test (norm(vec(x_ori)-vec(x_cc2))/norm(vec(x_ori))) < 2e-4

   #  Perform GCC from kspace
   kdataGCC,ccMat3 = geometricCoilCompression(kdata, 6)
   acqDataGCC = AcquisitionData(kdataGCC)
   x_gcc = reconstruction(acqDataGCC, params)
   x_gcc = mergeChannels(x_gcc)
   @test (norm(vec(x_gcc)-vec(x_ori))/norm(vec(x_ori))) < 3e-7


end

@testset "CoilCompression" begin
  test2DCC()
  test3DCC()
end
