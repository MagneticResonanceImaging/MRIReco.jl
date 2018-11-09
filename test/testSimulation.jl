# Basic Image(Shepp Logan Phantom) and

function test_kdata(N::Int64=32)
    @info "Testing simulating kdata with NFFT and exact evauluation"
    I = shepp_logan(N)
    tr = RadialTrajectory(N,N)
    println("Simulating kdata using NFFT")
    @time acqDataNFFT = simulation_fast(tr,I)
    @info "Simulating kdata rigorously"
    @time acqDataExplicit = simulation_explicit(tr,I)
    relError =  norm(acqDataExplicit.kdata[:]-acqDataNFFT.kdata[:]) / norm(acqDataExplicit.kdata[:])
    println("Relative error NFFT vs EXACT: ", relError)
    @test relError < 1e-2

end

function test_kdataMultipleSlices(N::Int64=32)
    @info "Testing simulating multiple 2d-slices with NFFT and exact evauluation"
    sh = ComplexF64.(shepp_logan(N))
    I = cat(sh,0.6*sh,0.3*sh,dims=3)
    tr = SimpleCartesianTrajectory(N,N)
    println("Simulating kdata using NFFT")
    @time acqDataNFFT = simulation(tr,I,opName="fast")
    @info "Simulating kdata rigorously"
    @time acqDataExplicit = simulation(tr,I,opName="explicit")
    relError =  norm(acqDataExplicit.kdata[:]-acqDataNFFT.kdata[:]) / norm(acqDataExplicit.kdata[:])
    println("Relative error NFFT vs EXACT: ", relError)
    @test relError < 1e-2

end

function test_kdata3d(N::Int64=32)
    @info "Testing simulating 3d-kdata with NFFT and exact evauluation"
    sh = ComplexF64.(shepp_logan(N))
    I = cat(sh,0.6*sh,0.3*sh,dims=3)
    tr = CartesianTrajectory3D(N,N,numSlices=3)
    println("Simulating kdata using NFFT")
    @time acqDataNFFT = simulation(tr,I,opName="fast")
    @info "Simulating kdata rigorously"
    @time acqDataExplicit = simulation(tr,I,opName="explicit")
    relError =  norm(acqDataExplicit.kdata[:]-acqDataNFFT.kdata[:]) / norm(acqDataExplicit.kdata[:])
    println("Relative error NFFT vs EXACT: ", relError)
    @test relError < 1e-2

end

function test_kdataWithCorrection(N::Int64=32)
    # Testing generating kdata with fieldinhomogeneity
    I = shepp_logan(N)
    @info "Testing simulating kdata with correctionterm"
    fmap = quadraticFieldmap(N,N,125*2pi)[:,:,1]
    rmap = relaxationMap(I,5.0,50.0)
    cmap = rmap + 1im*fmap
    tr = SpiralTrajectory(N,2*N,TE=0.0,AQ=32e-3)

    println("Simulating kdata using NFFT to approx. correctionterm ...")
    @time acqDataNFFT = simulation_fast(tr,I,cmap;method="nfft")
    println("Simulating kdata using Least Squares to approx. correctionterm...")
    @time acqDataLSQR = simulation_fast(tr,I,cmap;method="leastsquare")
    println("Simulating kdata rigorously...")
    @time acqDataExplicit = simulation_explicit(tr,I,cmap)
    # Calculating and testing the relative error of kdata
    relErrorLeastSquare = norm(acqDataExplicit.kdata[:]-acqDataLSQR.kdata[:]) / norm(acqDataExplicit.kdata[:])
    println("Relative Error of leastsquare method: ", relErrorLeastSquare)
    @test relErrorLeastSquare < 1e-3

    relErrorNFFT = norm(acqDataExplicit.kdata[:] - acqDataNFFT.kdata[:]) / norm(acqDataExplicit.kdata[:])
    println("Relative Error of leastsquare method: ", relErrorNFFT)
    @test relErrorNFFT < 1e-3

end


function testSimulation()
  @testset "simulations" begin
    test_kdata()
    test_kdataMultipleSlices()
    # test_kdata3d()
    test_kdataWithCorrection(16)
  end
end
