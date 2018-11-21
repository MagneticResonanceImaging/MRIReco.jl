# test gridding reco
function testGriddingReco(N=32)
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
  params[:reco] = "direct"
  params[:shape] = (N,N)

  x_approx = reconstruction(acqData, params)
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-2
end

function testGriddingReco3d(N=32)
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
  params[:reco] = "direct"
  params[:shape] = (N,N,8)

  x_approx = vec(reconstruction(acqData, params))
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-2
end

# test CS reco
function testCSReco(N=32,redFac=1.1)
  # image
  x = shepp_logan(N)

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N

  acqData = simulation(x, params)
  acqData = MRIReco.sample_kspace(acqData, redFac, "poisson", calsize=5)
  acqData = convertUndersampledData(acqData)

  # reco
  params[:reco] = "standard"    # encoding model
  params[:shape] = (N,N)
  params[:sparseTrafoName] = "Wavelet" #sparse trafo
  params[:regularization] = "L1"       # regularization
  params[:lambdL1] = 1.e-3
  params[:solver] = "admm"    # solver
  params[:iterations] = 1000
  params[:ρ] = 1.0e-1

  x_approx = reconstruction(acqData, params)
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-1
end

# test CSSense Reco
function testCSSenseReco(N=32,redFac=1.1)
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
  params[:numSamplingPerProfile] = N

  acqData = simulation( real(x.*reshape(sensMaps[:,1,1],N,N)), params )
  acqData2 = simulation( real(x.*reshape(sensMaps[:,2,1],N,N)), params )
  acqData.kdata = cat(acqData.kdata, acqData2.kdata, dims=1)
  acqData.numCoils = 2
  acqData.samplePointer = [1,N*N+1]
  acqData = MRIReco.sample_kspace(acqData, redFac, "poisson", calsize=5)
  acqData = convertUndersampledData(acqData)

  # reco
  params[:reco] = "multiCoil"
  params[:shape] = (N,N)
  params[:senseMaps] = sensMaps
  params[:sparseTrafoName] = "Wavelet" # sparse trafo
  params[:regularization] = "L1"       # regularization
  params[:lambdL1] = 1.e-3
  params[:solver] = "admm"
  params[:iterations] = 1000
  params[:ρ] = 1.0e-1

  x_approx = vec(reconstruction(acqData, params))
  @test (norm(vec(x)-x_approx)/norm(vec(x))) < 1e-1
end

function testOffresonanceReco(N = 128)

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
  params[:reco] = "standard"
  params[:regularization] = "L2"
  params[:iterations] = 3
  params[:solver] = "admm"
  params[:shape] = (N,N)
  params[:correctionMap] = cmap

  Ireco = reconstruction(acqData, params)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1e-1
end


function testSENSEReco(N = 64)
  numCoils = 8
  I = shepp_logan(N)
  I = circularShutterFreq!(I,1)

  coilsens = birdcageSensitivity(N, 8, 1.5)

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
  params[:reco] = "multiCoil" #"standard"
  params[:shape] = (N,N)
  params[:regularization] = "L2"
  params[:iterations] = 3
  params[:solver] = "admm"
  params[:senseMaps] = reshape(coilsens, N*N, numCoils, 1)

  Ireco = reconstruction(acqData, params)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1e-1
end


function testOffresonanceSENSEReco(N = 64)
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
  params[:reco] = "multiCoil" #"standard"
  params[:shape] = (N,N)
  params[:regularization] = "L2"
  params[:iterations] = 3
  params[:solver] = "admm"
  params[:senseMaps] = reshape(coilsens, N*N, numCoils, 1)
  params[:correctionMap] = cmap

  Ireco = reconstruction(acqData, params)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1.6e-1
end






function testReco(N=32)
  @testset "Reconstructions" begin
    testGriddingReco()
    testGriddingReco3d()
    testCSReco()
    testCSSenseReco()
    testOffresonanceReco()
    testSENSEReco()
    testOffresonanceSENSEReco()
  end
end
