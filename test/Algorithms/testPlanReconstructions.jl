# test gridding reco
function testPlanGriddingReco(N=32;arrayType = Array)
  # image
  x = shepp_logan(N)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N

  acqData = simulation(x, params)

  #reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("direct")
  setAll!(plan, params)
  x_approx = reconstruct(build(plan), acqData)
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-2
end

# test gridding reco
function testPlanConvertKspace(N=32;arrayType = Array)
    # image
    x = shepp_logan(N)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N

    acqData = simulation(x, params)

    #reco
    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N)
    params[:arrayType] = arrayType

    plan = MRIRecoPlan("direct")
    setAll!(plan, params)
    x_approx = reconstruct(build(plan), acqData)
    kspace = kDataCart(acqData)
    x_approx2 = ifftshift(ifft(ifftshift(kspace)))
    diff = abs.(x_approx2[:,:,1,1,1,1]) .- abs.(x)
    any(diff .< 1e-10)
    @test (norm(vec(x)-vec(abs.(x_approx2[:,:,1,1,1,1])))/norm(vec(x))) < 1e-10

    ## undersample kspace
    T = eltype(kspace)
    mask_idx = sample_vdpoisson((N,N),2.0)
    mask = zeros(T,N,N,1,1,1,1)
    mask[mask_idx] .= 1

    kspace_cs = copy(kspace)
    kspace_cs = kspace_cs .* mask

    acqCS = AcquisitionData(kspace_cs)
    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N)
    params[:arrayType] = arrayType

    plan = MRIRecoPlan("direct")
    setAll!(plan, params)
    x_cs = reconstruct(build(plan), acqCS)
    x_cs2 = ifftshift(ifft(ifftshift(kspace_cs)))

    @test (norm(vec(abs.(x_cs))-vec(abs.(x_cs2[:,:,1,1,1,1])))/norm(vec(abs.(x_cs2)))) < 1e-10
end

function testPlanConvertKspace3D(N=32;arrayType = Array)
    # image
    x = shepp_logan(N)
    x = repeat(x,1,1,N)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian3D"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N
    params[:numSlices] = N

    acqData = simulation(x, params)

    #reco
    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N,N)
    params[:arrayType] = arrayType

    plan = MRIRecoPlan("direct")
    setAll!(plan, params)
    x_approx = reconstruct(build(plan), acqData)
    kspace = kDataCart(acqData)
    x_approx2 = ifftshift(ifft(ifftshift(kspace)))
    diff = abs.(x_approx2[:,:,:,1,1,1]) .- abs.(x)
    any(diff .< 1e-10)
    @test (norm(vec(x)-vec(abs.(x_approx2[:,:,:,1,1,1])))/norm(vec(x))) < 1e-10

    ## undersample kspace
    T = eltype(kspace)
    mask_idx = sample_vdpoisson((N,N),2.0)
    mask = zeros(T,N,N)
    mask[mask_idx] .= 1
    mask = repeat(mask,1,1,N,1,1,1)

    kspace_cs = copy(kspace)
    kspace_cs = kspace_cs .* mask

    acqCS = AcquisitionData(kspace_cs)
    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N,N)

    plan = MRIRecoPlan("direct")
    setAll!(plan, params)
    x_cs = reconstruct(build(plan), acqCS)
    x_cs2 = ifftshift(ifft(ifftshift(kspace_cs)))

    @test (norm(vec(abs.(x_cs))-vec(abs.(x_cs2[:,:,:,1,1,1])))/norm(vec(abs.(x_cs2)))) < 1e-10
end

function testPlanGriddingReco3d(N=32;arrayType = Array)
  sh = ComplexF64.(shepp_logan(N))
  x = cat(sh,0.9*sh,0.8*sh,0.7*sh,0.6*sh,0.5*sh,0.4*sh,0.3*sh,dims=3)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian3D"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N
  params[:numSlices] = 8

  acqData = simulation( real(x), params )

  # reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N,8)
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("direct")
  setAll!(plan, params)
  x_approx = vec(reconstruct(build(plan), acqData))
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-2
end

# test CS reco
function testPlanCSReco(N=32,redFac=1.1;sampling="poisson",arrayType = Array)
  # image
  x = shepp_logan(N)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N

  acqData = simulation(x, params)
  Random.seed!(1234)
  acqData = sample_kspace(acqData, redFac, sampling, calsize=5, profiles=false)
  acqData = addNoise(acqData,25.0)

  # reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:sparseTrafo] = "Wavelet" #sparse trafo
  params[:reg] = L1Regularization(1.e-3)       # regularization
  params[:solver] = ADMM    # solver
  params[:iterations] = 1000
  params[:rho] = 1.0e-1
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("standard")
  setAll!(plan, params)
  x_approx = reconstruct(build(plan), acqData)
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-1
end

function testPlanCSRecoMultCoil(N=32;type = ComplexF64,arrayType = Array)
  # image
  x = shepp_logan(N)
  smaps = birdcageSensitivity(32,2,3.0)
  smaps[:,:,:,1] = 10*smaps[:,:,:,1]

  # convert to type
  x = convert.(type,x)
  smaps = convert.(type,smaps)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "SpiralVarDens"
  params[:numProfiles] = floor(Int64, N/4)
  params[:numSamplingPerProfile] = 2*N
  params[:senseMaps] = smaps

  acqData = simulation(x, params)
  @test(typeof(acqData.kdata[1,1,1][1])==type)

  # reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:sparseTrafo] = nothing #sparse trafo
  params[:reg] = TVRegularization(2.e-3, shape = (N,N))       # regularization
  params[:solver] = ADMM    # solver
  params[:iterations] = 100
  params[:rho] = 1.0e-1
  params[:absTol] = 1.e-5
  params[:relTol] = 1.e-4
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("standard")
  setAll!(plan, params)
  x_approx = reshape( reconstruct(build(plan), acqData), 32,32,2)

  @test(typeof(x_approx[1])==type)

  err1 = norm(vec(smaps[:,:,1,1]) .* vec(x)-vec(x_approx[:,:,1]) )/norm(vec(smaps[:,:,1,1]) .* vec(x))
  err2 = norm(vec(smaps[:,:,1,2]) .* vec(x)-vec(x_approx[:,:,2]) )/norm(vec(smaps[:,:,1,2]) .* vec(x))

  @test err1 < 1.1e-1
  @test err2 < 1.1e-1 # This error increased. Why?
end

# test CSSense Reco
function testPlanCSSenseReco(N=32,redFac=1.1;arrayType = Array)
  # image
  x = shepp_logan(N)

  # coil sensitivites
  sensMaps = zeros(ComplexF64,N*N,2,1)
  sensMaps[1:floor(Int64, N*N/2),1,1] .= 1.0
  sensMaps[floor(Int64, N*N/2)+1:end,2,1] .= 1.0

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:senseMaps] = reshape(sensMaps,N,N,1,2)
  params[:numSamplingPerProfile] = N

  acqData = simulation(x,params)
  Random.seed!(1234)
  acqData = sample_kspace(acqData, redFac, "poisson", calsize=5)

  # reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:senseMaps] = reshape(sensMaps,N,N,1,2)
  params[:sparseTrafo] = "Wavelet" # sparse trafo
  params[:reg] = L1Regularization(1.e-3)       # regularization
  params[:solver] = ADMM
  params[:iterations] = 1000
  params[:rho] = 1.0e-1
  params[:absTol] = 1.e-5
  params[:relTol] = 1.e-4
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("multiCoil")
  setAll!(plan, params)
  x_approx = vec(reconstruct(build(plan), acqData))
  @test (norm(vec(x)-x_approx)/norm(vec(x))) < 1e-1
end

function testPlanCSRecoMultiCoilCGNR(;N=32,redFac=1.1,type = ComplexF32,arrayType = Array)
  x = shepp_logan(N)

  # coil sensitivites
  sensMaps = zeros(ComplexF64,N*N,2,1)
  sensMaps[1:floor(Int64, N*N/2),1,1] .= 1.0
  sensMaps[floor(Int64, N*N/2)+1:end,2,1] .= 1.0

  # convert to type
  x = convert.(type,x)
  sensMaps = convert.(type,sensMaps)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:senseMaps] = reshape(sensMaps,N,N,1,2)
  params[:numSamplingPerProfile] = N

  acqData = simulation(x,params)
  Random.seed!(1234)
  acqData = sample_kspace(acqData, redFac, "poisson", calsize=5)

  # reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:senseMaps] = reshape(sensMaps,N,N,1,2)
  params[:reg] = L2Regularization(1.e-3)       # regularization
  params[:solver] = CGNR
  params[:iterations] = 1
  params[:rho] = 1.0e-1
  params[:absTol] = 1.e-5
  params[:relTol] = 1.e-4
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("multiCoil")
  setAll!(plan, params)
  x_approx = vec(reconstruct(build(plan), acqData))
  # This test is expected to have a higher error since CGNR with CS will not work
  # We still keep this test since it revealed the issue #110
  @test (norm(vec(x)-x_approx)/norm(vec(x))) < 2e-1
end

function testPlanOffresonanceReco(N = 128; accelMethod="nfft",arrayType = Array)

  I = shepp_logan(N)
  I = circularShutterFreq!(I,1)
  cmap = 1im*quadraticFieldmap(N,N,125*2pi)

  # simulation parameters
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Spiral"
  params[:numProfiles] = 1
  params[:numSamplingPerProfile] = N*N
  params[:windings] = div(N,2)
  params[:AQ] = 3.0e-2
  params[:correctionMap] = cmap[:,:,1]

  # do simulation
  acqData = simulation(I, params)

  # reco parameters
  params = Dict{Symbol, Any}()
  params[:reg] = L2Regularization(0.0)
  params[:iterations] = 3
  params[:solver] = ADMM
  params[:reconSize] = (N,N)
  params[:correctionMap] = cmap
  params[:method] = accelMethod
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("standard")
  setAll!(plan, params)
  Ireco = reconstruct(build(plan), acqData)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1e-1
end

function testPlanSENSEReco(N = 64, T = ComplexF64;arrayType = Array)
  numCoils = 8
  I = T.(shepp_logan(N))
  I = circularShutterFreq!(I,1)

  coilsens = T.(birdcageSensitivity(N, 8, 1.5))

  # simulation parameters
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Spiral"
  params[:numProfiles] = 1
  params[:numSamplingPerProfile] = div(N*N,2)
  params[:windings] = div(N,4)
  params[:AQ] = 3.0e-2
  params[:senseMaps] = coilsens

  # do simulation
  acqData = simulation(I, params)

  # reco parameters
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:reg] = L2Regularization(0.0)
  params[:iterations] = 50
  params[:solver] = CGNR
  params[:senseMaps] = coilsens
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("multiCoil")
  setAll!(plan, params)
  Ireco = reconstruct(build(plan), acqData)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1e-1
end

function testPlanSENSEnoiseUnCorr(N = 64, T = ComplexF64;arrayType = Array)
  numCoils = 8
  R = 4
  Img = T.(shepp_logan(N))
  coilsens = T.(birdcageSensitivity(N, numCoils, 1.5))

  # simulation parameters
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = div(N,R)
  params[:numSamplingPerProfile] = N
  params[:senseMaps] = coilsens

  # do simulation
  acqData = simulation(Img, params)

  # Generate correlated noise
  psi = Matrix{Complex{Float64}}(I, numCoils, numCoils)
  psi[2:numCoils,1] .= 1
  noise = randn(Complex{Float64}, (2*N, numCoils))
  noise = noise * psi

  # add noise correlation to data
  noise_data = randn(Complex{Float64}, (N, div(N,R), numCoils))
  noise_data = reshape(noise_data, (N .* div(N,R), numCoils)) * psi
  acqData.kdata[1,1,1] = acqData.kdata[1,1,1] .+ noise_data

  # reco parameters
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:reg] = L2Regularization(0.0)
  params[:iterations] = 150
  params[:solver] = CGNR
  params[:arrayType] = arrayType
  params[:senseMaps] = coilsens


  plan = MRIRecoPlan("multiCoil")
  setAll!(plan, params)
  Ireco = reconstruct(build(plan), acqData)

  params[:noiseData] = noise
  setAll!(plan, params)
  IrecoUnCorr = reconstruct(build(plan), acqData)

  @test (norm(vec(Img)-vec(Ireco))/norm(vec(Img)-vec(IrecoUnCorr))) > 1
end

function testPlanOffresonanceSENSEReco(N = 64;arrayType = Array)
  numCoils = 8
  I = shepp_logan(N)
  I = circularShutterFreq!(I,1)

  coilsens = birdcageSensitivity(N, 8, 1.5)
  cmap = 1im*quadraticFieldmap(N,N,125*2pi)

  # simulation parameters
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Spiral"
  params[:numProfiles] = 1
  params[:numSamplingPerProfile] = div(N*N,2)
  params[:windings] = div(N,4)
  params[:AQ] = 3.0e-2
  params[:senseMaps] = coilsens
  params[:correctionMap] = cmap[:,:,1]

  # do simulation
  acqData = simulation(I, params)

  # reco parameters
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:reg] = L2Regularization(0.0)
  params[:iterations] = 3
  params[:solver] = ADMM
  params[:senseMaps] = coilsens
  params[:correctionMap] = cmap
  params[:arrayType] = arrayType
  params[:senseMaps] = coilsens

  plan = MRIRecoPlan("multiCoil")
  setAll!(plan, params)
  Ireco = reconstruct(build(plan), acqData)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1.6e-1
end

function testPlanSenseMultiEcho(N=32,T = ComplexF32;arrayType = Array)
  # image
  x = T.(shepp_logan(N))
  rmap = 20.0*ones(N,N)

  coilsens = T.(birdcageSensitivity(N, 8, 1.5))

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N
  params[:r2map] = rmap
  params[:T_echo] = [2.e-2, 4.e-2]
  params[:seqName] = "ME"
  params[:refocusingAngles] = Float64[pi,pi]
  params[:senseMaps] = coilsens

  acqData = simulation( real(x), params )

  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:reg] = L2Regularization(1.e-3)
  params[:iterations] = 1
  params[:solver] = CGNR
  params[:senseMaps] = reshape(coilsens,N,N,1,8)
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("multiCoil")
  setAll!(plan, params)
  x_approx = reshape(reconstruct(build(plan),acqData),N,N,:)

  # Reco multiCoilMultiEcho
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:reg] = L2Regularization(1.e-3)
  params[:iterations] = 1
  params[:solver] = CGNR
  params[:senseMaps] = reshape(coilsens,N,N,1,8)

  plan = MRIRecoPlan("multiCoilMultiEcho")
  setAll!(plan, params)
  x_approx2= reshape(reconstruct(build(plan), acqData),N,N,:)

  relErrorEcho1 = norm(x_approx - x_approx2)/norm(x_approx)
  @test relErrorEcho1 < 1e-6
end

function testPlanSenseMultiEchoMeasCoils(N=32;T = ComplexF32,arrayType = Array)
    # image
    x = T.(shepp_logan(N))
    rmap = 20.0*ones(N,N)

    coilsens = T.(measured2DSensitivity(N, 8))

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N
    params[:r2map] = rmap
    params[:T_echo] = [2.e-2, 4.e-2]
    params[:seqName] = "ME"
    params[:refocusingAngles] = Float64[pi,pi]
    params[:senseMaps] = coilsens

    acqData = simulation( real(x), params )

    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N)
    params[:reg] = L2Regularization(1.e-3)
    params[:iterations] = 1
    params[:solver] = CGNR
    params[:senseMaps] = reshape(coilsens,N,N,1,8)
    params[:arrayType] = arrayType

    plan = MRIRecoPlan("multiCoil")
    setAll!(plan, params)
    x_approx = reshape(reconstruct(build(plan), acqData),N,N,:)

    # Reco multiCoilMultiEcho
    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N)
    params[:reg] = L2Regularization(1.e-3)
    params[:iterations] = 1
    params[:solver] = CGNR
    params[:senseMaps] = reshape(coilsens,N,N,1,8)
    params[:arrayType] = arrayType

    plan = MRIRecoPlan("multiCoilMultiEcho")
    setAll!(plan, params)
    x_approx2= reshape(reconstruct(build(plan), acqData),N,N,:)

    relErrorEcho1 = norm(x_approx - x_approx2)/norm(x_approx)
    @test relErrorEcho1 < 1e-6
end

function testPlanRecoMultiEcho(N=32;arrayType = Array)
  # image
  x = ComplexF64.(shepp_logan(N))
  rmap = 20.0*ones(N,N)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N
  params[:r2map] = rmap
  params[:T_echo] = [2.e-2, 4.e-2]
  params[:seqName] = "ME"
  params[:refocusingAngles] = Float64[pi,pi]

  acqData = simulation( real(x), params )

  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)

  plan = MRIRecoPlan("direct")
  setAll!(plan, params)
  x_approx = reshape(reconstruct(build(plan), acqData),N,N,2)

  relErrorEcho1 = norm(exp(-20.0*2.e-2)*x - x_approx[:,:,1])/norm(exp(-20.0*2.e-2)*x)
  @test relErrorEcho1 < 1e-3
  relErrorEcho2 = norm(exp(-20.0*4.e-2)*x - x_approx[:,:,2])/norm(exp(-20.0*4.e-2)*x)
  @test relErrorEcho2 < 1e-3

  # Reco multiEcho
  params = Dict{Symbol, Any}()
  params[:reconSize] = (N,N)
  params[:reg] = L2Regularization(1.e-3)
  params[:iterations] = 1
  params[:solver] = CGNR
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("multiEcho")
  setAll!(plan, params)
  x_approx2 = reshape(reconstruct(build(plan), acqData),N,N,:)

  relErrorEcho3 = norm(x_approx - x_approx2)/norm(x_approx)
  @test relErrorEcho3 < 1e-3
end

function testPlanCSReco3d(N=64;arrayType = Array)
  sh = ComplexF64.(shepp_logan(N))
  I = cat(sh,0.9*sh,0.8*sh,0.7*sh,0.6*sh,0.5*sh,0.4*sh,0.3*sh,dims=3)
  I = permutedims(I,[3,1,2])

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian3D"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = 8
  params[:numSlices] = N

  acqData = simulation( real(I), params )
  Random.seed!(1234)
  acqData = sample_kspace(acqData,1.5,"poisson",calsize=15)

  # 3d reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (8,N,N)
  params[:sparseTrafo] = "Wavelet" #"nothing" #sparse trafo
  params[:reg] = L1Regularization(1.e-3) #"TV"       # regularization
  params[:solver] = ADMM    # solver
  params[:iterations] = 1000
  params[:rho] = 1.0e-1
  params[:absTol] = 1.e-4
  params[:relTol] = 1.e-4
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("standard")
  setAll!(plan, params)
  Ireco = reconstruct(build(plan), acqData)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 2e-1
end

function testPlanCSSenseReco3d(N=128;arrayType = Array)
  sh = ComplexF64.(shepp_logan(N))
  I = cat(sh,0.9*sh,0.8*sh,0.7*sh,0.6*sh,0.5*sh,0.4*sh,0.3*sh,dims=3)
  I = permutedims(I,[3,1,2])

  sensMaps = zeros(ComplexF64,8,N*N,2)
  sensMaps[:,1:floor(Int64, N*N/2),1] .= 1.0
  sensMaps[:,floor(Int64, N*N/2)+1:end,2] .= 1.0

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian3D"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = 8
  params[:numSlices] = N
  params[:senseMaps] = reshape(sensMaps,8,N,N,2)

  acqData = simulation( real(I), params )
  Random.seed!(1234)
  acqData = sample_kspace(acqData,4.0,"poisson",calsize=15)

  # 3d reco
  params = Dict{Symbol, Any}()
  params[:reconSize] = (8,N,N)
  params[:senseMaps] = reshape(sensMaps,8,N,N,2)
  params[:sparseTrafo] = "Wavelet" #"nothing" #sparse trafo
  params[:reg] = L1Regularization(1.e-3) #"TV"       # regularization
  params[:solver] = ADMM    # solver
  params[:iterations] = 1000
  params[:rho] = 1.0e-1
  params[:absTol] = 1.e-5
  params[:relTol] = 1.e-4
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("multiCoil")
  setAll!(plan, params)
  Ireco = reconstruct(build(plan), acqData)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1e-1
end

function testPlanRegridding(N=64;arrayType = Array)

  # image
  x = shepp_logan(N)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Radial"
  params[:numProfiles] = round(Int64,1.25*N)
  params[:numSamplingPerProfile] = round(Int64,1.25*N)

  acqDataRad = simulation(x, params)

  params[:trajName] = "Cartesian"
  params[:numProfiles] = N
  params[:numSamplingPerProfile] = N
  acqDataCart = simulation(x, params)

  # regridding
  params = Dict{Symbol, Any}()
  acqDataReg = regrid(acqDataRad,(N,N))
  # cartesian reconstruction
  params[:reconSize] = (N,N)
  params[:solver] = CGNR
  params[:reg] = L2Regularization(0.0)
  params[:iterations] = 3
  params[:arrayType] = arrayType

  plan = MRIRecoPlan("standard")
  setAll!(plan, params)
  x_reg = collect(reshape(reconstruct(build(plan), acqDataReg),N,N))
  circularShutter!(x_reg)
  x_rad = collect(reshape(reconstruct(build(plan), acqDataRad),N,N))
  circularShutter!(x_rad)

  k_reg = reshape(acqDataReg.kdata[1],N,N)
  k_cart = reshape(acqDataCart.kdata[1],N,N)
  circularShutter!(k_reg)
  circularShutter!(k_cart)

  @test (norm(vec(k_cart)-vec(k_reg))/norm(vec(k_cart))) < 1e-1
  @test (norm(vec(x_rad).-vec(x_reg))/norm(vec(x_rad))) < 1e-2
end

function testPlanRecos(N=32;arrayType = Array)
    @testset "Reconstructions: $arrayType" begin
        @testset "MultiEcho" begin
            @testset testPlanRecoMultiEcho(;arrayType)
            @testset testPlanSenseMultiEcho(N;arrayType)
            @testset testPlanSenseMultiEchoMeasCoils(N;arrayType)
        end

        @testset "Convert k-space" begin
            @testset "K-Space" testPlanConvertKspace(;arrayType)
            @testset "K-Space 3D" testPlanConvertKspace3D(;arrayType)
        end

        @testset "Gridding" begin
            @testset "Gridding" testPlanGriddingReco(;arrayType)
            @testset "Gridding 3D" testPlanGriddingReco3d(;arrayType)
            @testset "Regridding" testPlanRegridding(;arrayType)
        end

        @testset "CS" begin
            sampling = ["random", "poisson", "vdPoisson"] # "lines"
            @testset "Sampling" for samp in sampling
                testPlanCSReco(sampling=samp;arrayType)
            end
            @testset "MultiCoil F64" testPlanCSRecoMultCoil(type = ComplexF64;arrayType)
            @testset "MultiCoil F32" testPlanCSRecoMultCoil(type = ComplexF32;arrayType)
            @testset "Sense" testPlanCSSenseReco(;arrayType)
            @testset "3D" testPlanCSReco3d(;arrayType)
            @testset "Sense 3D" testPlanCSSenseReco3d(;arrayType)
            @testset "Multi Coil CGNR" testPlanCSRecoMultiCoilCGNR(type = ComplexF64;arrayType)
            @testset "Multi Coil CGNR" testPlanCSRecoMultiCoilCGNR(type = ComplexF32;arrayType)
        end

        @testset "off-resonance" begin
            accelMethods = ["nfft", "hist", "leastsquare"]
            for a in accelMethods
                !Sys.iswindows() && testPlanOffresonanceReco(accelMethod=a;arrayType)
                !Sys.iswindows() && testPlanOffresonanceSENSEReco(;arrayType)
            end
        end

        @testset "SENSE" begin
            @testset "SENSE F32" testPlanSENSEReco(64, ComplexF32;arrayType)
            @testset "SENSE F64" testPlanSENSEReco(64, ComplexF64;arrayType)
            @testset "SENSE F64 Uncorr" testPlanSENSEnoiseUnCorr(64, ComplexF64;arrayType)
        end
    end
end

for arrayType in arrayTypes 
  testPlanRecos(;arrayType)
end