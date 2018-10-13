# Basic Image(Shepp Logan Phantom) and

function test_kdata(N::Int64=32)
    @info "Testing simulating kdata with NFFT and exact evauluation"
    I = shepp_logan(N)
    tr = RadialTrajectory(N,N)
    println("Simulating kdata using NFFT")
    @time aqDataNFFT = simulation_nfft(tr,I)
    @info "Simulating kdata rigorously"
    @time aqDataExplicit = simulation_explicit(tr,I)
    relError =  norm(aqDataExplicit.kdata[:]-aqDataNFFT.kdata[:]) / norm(aqDataExplicit.kdata[:])
    println("Relative error NFFT vs EXACT: ", relError)
    @test relError < 1e-2

end

function test_kdataWithCorrection(N::Int64=32)
    # Testing generating kdata with fieldinhomogeneity
    I = shepp_logan(N)
    @info "Testing simulating kdata with correctionterm"
    fmap = quadraticFieldmap(N,N,125*2pi)
    rmap = relaxationMap(I,5.0,50.0)
    cmap = rmap + 1im*fmap
    tr = SpiralTrajectory(N,2*N,TE=0.0,AQ=32e-3)

    println("Simulating kdata using NFFT to approx. correctionterm ...")
    @time aqDataNFFT = simulation_nfft(tr,I,cmap;method="nfft")
    println("Simulating kdata using Least Squares to approx. correctionterm...")
    @time aqDataLSQR = simulation_nfft(tr,I,cmap;method="leastsquare")
    println("Simulating kdata rigorously...")
    @time aqDataExplicit = simulation_explicit(tr,I,cmap)
    # Calculating and testing the relative error of kdata
    relErrorLeastSquare = norm(aqDataExplicit.kdata[:]-aqDataLSQR.kdata[:]) / norm(aqDataExplicit.kdata[:])
    println("Relative Error of leastsquare method: ", relErrorLeastSquare)
    @test relErrorLeastSquare < 1e-3

    relErrorNFFT = norm(aqDataExplicit.kdata[:] - aqDataNFFT.kdata[:]) / norm(aqDataExplicit.kdata[:])
    println("Relative Error of leastsquare method: ", relErrorNFFT)
    @test relErrorNFFT < 1e-3

end


function testSimulation()
  @testset "simulations" begin
    test_kdata()
    test_kdataWithCorrection(16)
  end
end
