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

  acqDataCC, smapsCC = geometricCC_2d(acqData, smaps, 6)
  x_cc = reconstruction(acqDataCC, params)
  x_cc = mergeChannels(x_cc)

  @test (norm(vec(x_ori)-vec(x_cc))/norm(vec(x_ori))) < 1e-1
end

@testset "CoilCompression" begin
    test2DCC()
end
