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

function test2DMSCC(N=64,nSl = 3)
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

  acqData_MS = deepcopy(acqData)
  acqData_MS.kdata = [acqData.kdata[1,1,1] for echo=1:1, sl=1:nSl, rep=1:1]

  #reconstruction
  params = Dict{Symbol,Any}()
  params[:reco] = "direct"
  params[:reconSize] = acqData_MS.encodingSize

  x_ori = reconstruction(acqData_MS, params)
  x_ori = mergeChannels(x_ori)

  acqDataCC,ccMat = softwareCoilCompression(acqData_MS, 6)
  x_cc = reconstruction(acqDataCC, params)
  x_cc = mergeChannels(x_cc)
  @test (norm(vec(x_ori)-vec(x_cc))/norm(vec(x_ori))) < 2e-4

  acqDataCC,ccMat = geometricCoilCompression(acqData_MS, 6)
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

   #  Perform GCC from acq
   acqDataGCC,ccMat3 = geometricCoilCompression(acqData, 6)
   x_gcc = reconstruction(acqDataGCC, params)
   x_gcc = mergeChannels(x_gcc)
   @test (norm(vec(x_gcc)-vec(x_ori))/norm(vec(x_ori))) < 3e-7
end

function simulate2D(img, smaps; nSl = 1, nContr = 1, nRep = 1)
  N = size(img, 1)
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = N
  params[:numSamplingPerProfile] = N
  params[:senseMaps] = smaps
  acqData = simulation(img, params)
  # replicate the data along contrasts, slices and repetitions with different scalings,
  # so that each block is distinct
  acqData.kdata = [acqData.kdata[1,1,1] .* (contr + 10sl + 100rep)
                   for contr = 1:nContr, sl = 1:nSl, rep = 1:nRep]
  acqData.traj = fill(acqData.traj[1], nContr)
  acqData.subsampleIndices = fill(acqData.subsampleIndices[1], nContr)
  return acqData
end

function testApplyCC(N = 32)
  img = shepp_logan(N)
  smaps = birdcageSensitivity(N, 8, 1.5)
  numVC = 4

  # acqData with several contrasts, slices and repetitions
  acqData = simulate2D(img, smaps; nSl = 3, nContr = 2, nRep = 2)
  acqDataCC, ccMat_vec = softwareCoilCompression(acqData, numVC; sContr = 2, sRep = 2)
  @test length(ccMat_vec) == 3
  @test size(ccMat_vec[1]) == (8, numVC)
  @test numChannels(acqDataCC) == numVC
  for contr = 1:2, sl = 1:3, rep = 1:2
    @test acqDataCC.kdata[contr, sl, rep] ≈ acqData.kdata[contr, sl, rep] * ccMat_vec[sl]
  end
  # the input is not modified and the other fields are kept
  @test numChannels(acqData) == 8
  @test acqDataCC.subsampleIndices == acqData.subsampleIndices
  @test acqDataCC.encodingSize == acqData.encodingSize

  # applyCoilCompression with a single matrix uses it for every slice
  acqDataCC1 = applyCoilCompression(acqData, ccMat_vec[1])
  @test acqDataCC1.kdata[1, 3, 2] ≈ acqData.kdata[1, 3, 2] * ccMat_vec[1]
  @test_throws DimensionMismatch applyCoilCompression(acqData, ccMat_vec[1:2])

  # kDataCart and applyCoilCompression commute
  kspace = kDataCart(acqData)
  @test applyCoilCompression(kspace, ccMat_vec) ≈ kDataCart(acqDataCC)

  # softwareCoilCompression on the k-space compresses all contrasts and repetitions
  kspaceCC, ccMat = softwareCoilCompression(kspace, numVC; sContr = 2, sRep = 2)
  @test size(kspaceCC) == (N, N, 3, numVC, 2, 2)
  @test kspaceCC ≈ applyCoilCompression(kspace, ccMat)
  @test kspaceCC[:, :, :, :, 1, 1] ≈ reshape(reshape(kspace[:, :, :, :, 1, 1], :, 8) * ccMat, N, N, 3, numVC)
  @test_throws DimensionMismatch applyCoilCompression(kspace, ccMat_vec[1:2])

  # Compressed sensitivity maps are the sensitivities of the virtual coils:
  # simulating with them gives the compressed data
  acqData1 = simulate2D(img, smaps)
  acqData1CC, ccMat1 = softwareCoilCompression(acqData1, numVC)
  smapsCC = applyCoilCompressionSensitivityMaps(smaps, ccMat1[1])
  @test size(smapsCC) == (N, N, 1, numVC)
  @test applyCoilCompressionSensitivityMaps(smaps, ccMat1) ≈ smapsCC
  acqDataSimCC = simulate2D(img, smapsCC)
  @test acqDataSimCC.kdata[1,1,1] ≈ acqData1CC.kdata[1,1,1]

  # one matrix per slice for multi-slice maps
  smapsMS = cat(smaps, 2 .* smaps, 3 .* smaps, dims = 3)
  smapsMSCC = applyCoilCompressionSensitivityMaps(smapsMS, ccMat_vec)
  for sl = 1:3
    @test vec(smapsMSCC[:, :, sl, :]) ≈ vec((sl .* reshape(smaps, N*N, 8)) * ccMat_vec[sl])
  end
  @test_throws DimensionMismatch applyCoilCompressionSensitivityMaps(smapsMS, ccMat_vec[1:2])
end

function testGCCInvalidDim()
  kspace = zeros(ComplexF32, 8, 8, 1, 4, 1, 1)
  for dim in (-1, 0, 4)
    @test_throws ArgumentError geometricCoilCompression(kspace, 2; dim = dim)
  end
end

# Reference implementation of GCC: compress each position along `dim` of the hybrid space explicitly
centeredIFFT(x, dim) = fftshift(ifft(ifftshift(x, dim), dim), dim)
centeredFFT(x, dim) = fftshift(fft(ifftshift(x, dim), dim), dim)
function gccReference(kspace, ccMat, dim)
  hybrid = centeredIFFT(kspace, dim)
  nChan, numVC = size(ccMat, 2), size(ccMat, 3)
  out = zeros(eltype(kspace), size(kspace)[1:3]..., numVC, size(kspace)[5:6]...)
  for p in 1:size(kspace, dim)
    h = collect(selectdim(hybrid, dim, p:p))
    sz = size(h)
    c = reshape(permutedims(reshape(h, prod(sz[1:3]), nChan, prod(sz[5:6])), (1, 3, 2)), :, nChan) * ccMat[p, :, :]
    selectdim(out, dim, p:p) .= reshape(permutedims(reshape(c, prod(sz[1:3]), prod(sz[5:6]), numVC), (1, 3, 2)),
                                        sz[1:3]..., numVC, sz[5:6]...)
  end
  return centeredFFT(out, dim)
end

function testGCCOddSizes()
  # odd sizes everywhere (issue 5): 7 channels, 3 virtual coils, 3 contrasts, odd matrix size
  kspace = randn(ComplexF64, 9, 7, 5, 7, 3, 1)
  for dim in 1:3
    kspaceCC, ccMat = geometricCoilCompression(kspace, 3; dim = dim, sContr = 2)
    @test size(ccMat) == (size(kspace, dim), 7, 3)
    @test size(kspaceCC) == (9, 7, 5, 3, 3, 1)
    @test kspaceCC ≈ gccReference(kspace, ccMat, dim)
    @test kspaceCC ≈ applyCoilCompression(kspace, ccMat; dim = dim)
    # each matrix has orthonormal columns
    @test ccMat[2, :, :]' * ccMat[2, :, :] ≈ I
    # the calibration uses contrast sContr = 2
    _, ccMat2 = geometricCoilCompression(kspace[:, :, :, :, 2:2, :], 3; dim = dim)
    @test ccMat ≈ ccMat2
  end
end

function testGCCDims()
  # compressing along dim 2 or 3 (issue 6) is the same as permuting that dimension to the front
  kspace = randn(ComplexF64, 6, 8, 10, 4, 1, 1)
  for (dim, perm) in ((2, (2, 1, 3, 4, 5, 6)), (3, (3, 1, 2, 4, 5, 6)))
    kspaceCC, ccMat = geometricCoilCompression(kspace, 2; dim = dim)
    kspaceCCp, ccMatp = geometricCoilCompression(permutedims(kspace, perm), 2; dim = 1)
    @test ccMat ≈ ccMatp
    @test kspaceCC ≈ permutedims(kspaceCCp, invperm(perm))
  end
end

function testGCCDefaults()
  # default number of virtual coils (issue 8): all channels, the data are unchanged
  kspace = randn(ComplexF64, 8, 8, 1, 4, 1, 1)
  kspaceCC, ccMat = geometricCoilCompression(kspace)
  @test size(ccMat) == (8, 4, 4)
  @test norm(kspaceCC) ≈ norm(kspace)
  @test_throws ArgumentError geometricCoilCompression(kspace, 5)
  @test_throws ArgumentError softwareCoilCompression(kspace, 0)
  # deprecated positional sContr / sRep
  @test_deprecated softwareCoilCompression(kspace, 2, 1, 1)
end

function testGCCAcqData(N = 32)
  img = shepp_logan(N)
  smaps = birdcageSensitivity(N, 8, 1.5)
  acqData = simulate2D(img, smaps; nSl = 2, nContr = 2)
  acqDataSub = sample_kspace(acqData, 2.0, "poisson", calsize = 8)

  # default number of virtual coils (issue 8)
  @test numChannels(first(geometricCoilCompression(acqData))) == 8

  # undersampled data stay undersampled (issue 9): the profiles are missing along dimension 2,
  # dimension 1 (readout) is fully sampled
  @test_throws ArgumentError geometricCoilCompression(acqDataSub, 4; dim = 2)
  for (dim, acq) in ((1, acqDataSub), (2, acqData))
    acqDataCC, ccMat = geometricCoilCompression(acq, 4; dim = dim)
    @test numChannels(acqDataCC) == 4
    @test acqDataCC.subsampleIndices == acq.subsampleIndices
    @test size(acqDataCC.kdata[1, 1, 1], 1) == size(acq.kdata[1, 1, 1], 1)
    @test kDataCart(acqDataCC) ≈ applyCoilCompression(kDataCart(acq), ccMat; dim = dim)
    @test kDataCart(acqDataCC) ≈ gccReference(kDataCart(acq), ccMat, dim)
  end
  @test_throws ArgumentError geometricCoilCompression(acqData, 4; dim = 3)

  # compressed sensitivity maps are the sensitivities of the virtual coils (issue 3)
  acqData1 = simulate2D(img, smaps)
  for dim in (1, 2)
    acqData1CC, ccMat = geometricCoilCompression(acqData1, 4; dim = dim)
    smapsCC = applyCoilCompressionSensitivityMaps(smaps, ccMat; dim = dim)
    @test size(smapsCC) == (N, N, 1, 4)
    @test simulate2D(img, smapsCC).kdata[1, 1, 1] ≈ acqData1CC.kdata[1, 1, 1]
  end
end

@testset "CoilCompression" begin
  testGCCInvalidDim()
  testGCCOddSizes()
  testGCCDims()
  testGCCDefaults()
  testGCCAcqData()
  test2DCC()
  test3DCC()
  test2DMSCC()
  testApplyCC()
end
