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
function testCSReco(N=32,redFac=1.1;sampling="poisson")
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
  acqData = MRIReco.sample_kspace(acqData, redFac, sampling, calsize=5)
  acqData = MRIReco.addNoise(acqData,25.0)

  # reco
  params[:reco] = "standard"    # encoding model
  params[:shape] = (N,N)
  params[:sparseTrafoName] = "Wavelet" #sparse trafo
  params[:regularization] = "L1"       # regularization
  params[:λ] = 1.e-3
  params[:solver] = "admm"    # solver
  params[:iterations] = 1000
  params[:ρ] = 1.0e-1

  x_approx = reconstruction(acqData, params)
  @test (norm(vec(x)-vec(x_approx))/norm(vec(x))) < 1e-1
end

function testCSRecoMultCoil(N=32)
  # image
  x = shepp_logan(N)
  smaps = birdcageSensitivity(32,2,3.0)
  smaps[:,:,:,1] = 10*smaps[:,:,:,1]

  # simulation
  params = Dict{Symbol, Any}()
  params[:simulation] = "fast"
  params[:trajName] = "SpiralVarDens"
  params[:numProfiles] = floor(Int64, N/4)
  params[:numSamplingPerProfile] = 2*N
  params[:senseMaps] = smaps

  acqData = simulation(x, params)

  # reco
  params[:reco] = "standard"    # encoding model
  params[:shape] = (N,N)
  params[:sparseTrafoName] = "nothing" #sparse trafo
  params[:regularization] = "TV"       # regularization
  params[:λ] = 2.e-3
  params[:solver] = "admm"    # solver
  params[:iterations] = 100
  params[:ρ] = 1.0e-1
  params[:absTol] = 1.e-5
  params[:relTol] = 1.e-4

  x_approx = reshape( reconstruction(acqData, params), 32,32,2)

  err1 = norm(vec(smaps[:,:,1,1]) .* vec(x)-vec(x_approx[:,:,1]) )/norm(vec(smaps[:,:,1,1]) .* vec(x))
  err2 = norm(vec(smaps[:,:,1,2]) .* vec(x)-vec(x_approx[:,:,2]) )/norm(vec(smaps[:,:,1,2]) .* vec(x))

  @test err1 < 1e-1
  @test err2 < 1e-1
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
  params[:senseMaps] = reshape(sensMaps,N,N,1,2)
  params[:numSamplingPerProfile] = N

  acqData = simulation(x,params)
  Random.seed!(1234)
  acqData = MRIReco.sample_kspace(acqData, redFac, "poisson", calsize=5)

  # reco
  params[:reco] = "multiCoil"
  params[:shape] = (N,N)
  params[:senseMaps] = sensMaps
  params[:sparseTrafoName] = "Wavelet" # sparse trafo
  params[:regularization] = "L1"       # regularization
  params[:λ] = 1.e-3
  params[:solver] = "admm"
  params[:iterations] = 1000
  params[:ρ] = 1.0e-1
  params[:absTol] = 1.e-5
  params[:relTol] = 1.e-4

  x_approx = vec(reconstruction(acqData, params))
  @test (norm(vec(x)-x_approx)/norm(vec(x))) < 1e-1
end

function testOffresonanceReco(N = 128; accelMethod="nfft")

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
  params[:method] = accelMethod

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

function testDirectRecoMultiEcho(N=32)
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

  params[:reco] = "direct"
  params[:shape] = (N,N)

  x_approx = reshape(reconstruction(acqData,params),N,N,2)

  relErrorEcho1 = norm(exp(-20.0*2.e-2)*x - x_approx[:,:,1])/norm(exp(-20.0*2.e-2)*x)
  @test relErrorEcho1 < 1e-3
  relErrorEcho2 = norm(exp(-20.0*4.e-2)*x - x_approx[:,:,2])/norm(exp(-20.0*4.e-2)*x)
  @test relErrorEcho2 < 1e-3
end

function testCSReco3d(N=64)
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
  acqData = MRIReco.sample_kspace(acqData,1.5,"poisson",calsize=5)

  # 3d reco
  params[:reco] = "standard"    # encoding model
  params[:shape] = (8,N,N)
  params[:sparseTrafoName] = "nothing" #sparse trafo
  params[:regularization] = "TV"       # regularization
  params[:λ] = 1.e-3
  params[:solver] = "admm"    # solver
  params[:iterations] = 1000
  params[:ρ] = 1.0e-1
  params[:absTol] = 1.e-4
  params[:relTol] = 1.e-3

  Ireco = reconstruction(acqData, params)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 1e-1
end

function testCSSenseReco3d(N=64)
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
  acqData = MRIReco.sample_kspace(acqData,4.0,"poisson",calsize=5)

  # 3d reco
  params[:reco] = "multiCoil"    # encoding model
  params[:shape] = (8,N,N)
  params[:senseMaps] = reshape(sensMaps,8,N,N,2)
  params[:sparseTrafoName] = "nothing" #sparse trafo
  params[:regularization] = "TV"       # regularization
  params[:λ] = 1.e-3
  params[:solver] = "admm"    # solver
  params[:iterations] = 1000
  params[:ρ] = 1.0e-1
  params[:absTol] = 1.e-5
  params[:relTol] = 1.e-4

  Ireco = reconstruction(acqData, params)

  @test (norm(vec(I)-vec(Ireco))/norm(vec(I))) < 5e-2
end

function testReco(N=32)
  @testset "Reconstructions" begin
    testGriddingReco()
    testGriddingReco3d()
    sampling = ["random", "poisson", "vdPoisson"] # "lines"
    for samp in sampling
      testCSReco(sampling=samp)
    end
    testCSRecoMultCoil()
    testCSSenseReco()
    testCSReco3d()
    testCSSenseReco3d()
    accelMethods = ["nfft", "hist", "leastsquare"]
    for a in accelMethods
      testOffresonanceReco(accelMethod=a)
    end
    testSENSEReco()
    testOffresonanceSENSEReco()
    testDirectRecoMultiEcho()
  end
end

testReco()
