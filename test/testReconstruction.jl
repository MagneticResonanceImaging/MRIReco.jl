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

  aqData = simulation(x, params)

  #reco
  params[:reco] = "nfft"
  params[:shape] = (N,N)

  x_approx = reconstruction(aqData, params)
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

  aqData = simulation(x, params)
  aqData = MRIReco.sample_kspace(aqData, redFac, "poisson", calsize=5)
  aqData = convertUndersampledData(aqData)

  # reco
  params[:reco] = "simple"    # encoding model
  params[:shape] = (N,N)
  params[:sparseTrafoName] = "Wavelet" #sparse trafo
  params[:regularization] = "L1"       # regularization
  params[:lambdL1] = 1.e-3
  params[:solver] = "admm"    # solver
  params[:iterations] = 1000
  params[:ρ] = 1.0e-1

  x_approx = reconstruction(aqData, params)
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-1
end

# test CSSense Reco
function testCSSenseReco(N=32,redFac=1.1)
  # image
  x = shepp_logan(N)

  # coil sensitivites
  sensMaps = zeros(N*N,2,1)
  sensMaps[1:floor(Int64, N*N/2),1,1] .= 1.0
  sensMaps[floor(Int64, N*N/2)+1:end,2,1] .= 1.0

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "Cartesian"
  params[:numProfiles] = floor(Int64, N)
  params[:numSamplingPerProfile] = N

  aqData = simulation( real(x.*reshape(sensMaps[:,1],N,N)), params )
  aqData2 = simulation( real(x.*reshape(sensMaps[:,2],N,N)), params )
  aqData.kdata = cat(aqData.kdata, aqData2.kdata, dims=1)
  aqData.numCoils = 2
  aqData.samplePointer = [1,N*N+1]
  aqData = MRIReco.sample_kspace(aqData, redFac, "poisson", calsize=5)
  aqData = convertUndersampledData(aqData)

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

  x_approx = vec(reconstruction(aqData, params))
  @test (norm(vec(x)-x_approx)/norm(vec(x))) < 1e-1
end


# test CSSense Reco
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
  params[:correctionMap] = cmap

  # do simulation
  aqData = simulation(I, params)

  # reco parameters
  params = Dict{Symbol, Any}()
  params[:reco] = "simple"
  params[:regularization] = "L2"
  params[:iterations] = 3
  params[:solver] = "admm"
  params[:shape] = (N,N)
  params[:cmap] = cmap
  params[:alpha] = 1.75
  params[:m] = 4.0
  params[:K] = 28
  Ireco = reconstruction(aqData, params)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1e-1
end

function testReco(N=32)
  @testset "Reconstructions" begin
    testGriddingReco()
    testCSReco()
    testCSSenseReco()
    testOffresonanceReco()
  end
end
