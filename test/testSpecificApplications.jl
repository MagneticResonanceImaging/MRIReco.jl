function testConversionAcqToRaw(N=32)
    # image
    x = shepp_logan(N)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = round(Int64,1.25*N)
    params[:numSamplingPerProfile] = round(Int64,1.25N)

    acqDataRad = simulation(x, params)

    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    params[:reconSize] = (N,N)

    Ireco = reconstruction(acqDataRad, params)

    # convert acqDataRad to rawRad
    rawRad = RawAcquisitionData(acqDataRad)
    acqRad = AcquisitionData(rawRad)

    Ireco2 = reconstruction(acqRad, params)

    @test (norm(vec(Ireco)-vec(Ireco2))/norm(vec(Ireco))) < 1e-5

    ## try that in 3D
    x = shepp_logan(N)
    x = repeat(x,1,1,N)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Kooshball"
    params[:numProfiles] = round(Int64,pi*N^2)
    params[:numSamplingPerProfile] = round(Int64,N)

    acqDataRad = simulation(x, params)

    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    params[:reconSize] = (N,N,N)

    Ireco = reconstruction(acqDataRad, params)

    # convert acqDataRad to rawRad
    rawRad = RawAcquisitionData(acqDataRad)
    acqRad = AcquisitionData(rawRad)

    Ireco2 = reconstruction(acqRad, params)
    @test (norm(vec(Ireco)-vec(Ireco2))/norm(vec(Ireco))) < 1e-4
end

function testBinning(N=32)
    # image
    x = shepp_logan(N)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Radial"
    params[:numProfiles] = round(Int64,900)
    params[:numSamplingPerProfile] = round(Int64,N)

    acqRad = simulation(x, params)

    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    params[:reconSize] = (N,N)

    # convert acqDataRad to rawRad
    rawRad = RawAcquisitionData(acqRad)

    # split into one acquisition with 6 times more projection
    rawRad2 = deepcopy(rawRad)
    for i in 1:length(rawRad.profiles)
        if  mod(i,3) == 0
            rawRad2.profiles[i].head.idx.contrast = 0
        else
            rawRad2.profiles[i].head.idx.contrast = 1
        end
    end
    ## edit metadata for trajectory
    rawRad2.params["trajectory"] = "custom"

    # transform to acq
    acqRad2 = AcquisitionData(rawRad2)

    @test length(acqRad2.kdata) == 2
    @test size(acqRad2.kdata[1],1)/N == 900/3
    @test size(acqRad2.kdata[2],1)/N == 2*900/3
    Ireco3 = reconstruction(acqRad2, params)

    @test (norm(vec(Ireco3[:,:,1,1])-vec(Ireco3[:,:,1,2]))/norm(vec(Ireco3[:,:,1,2]))) < 0.05
end

function testSubspace(N=32)
    T = ComplexF32
    x = T.(shepp_logan(N))
    rmap = 5.0*abs.(x)
    coilsens = T.(birdcageSensitivity(N, 8, 1.5))
    TEnum = collect(2.e-2:2.e-2:50.e-2)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N
    params[:r2map] = rmap
    params[:T_echo] = TEnum
    params[:seqName] = "ME"
    params[:refocusingAngles] = Float64.(repeat([pi],length(TEnum)))
    params[:senseMaps] = coilsens

    acqData = simulation( real(x), params )

    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N)
    params[:reco] = "multiCoilMultiEcho"
    #params[:sparseTrafo] = "nothing" #"nothing" #sparse trafo
    params[:reg] = L2Regularization(1.e-3)
    params[:iterations] = 1
    params[:solver] = CGNR
    params[:senseMaps] = reshape(coilsens,N,N,1,8)

    x_approx= reconstruction(acqData,params)

    ## create basis
    function createExpBasis(TE_vec::AbstractVector{T},T2_vec::AbstractVector{T}) where {T<:Real}
    nTE = length(TE_vec)
    nSimu = length(T2_vec)
    expSignal = zeros(T,nSimu,nTE)
    
    for (i,T2) in enumerate(T2_vec)
        expSignal[i,:] = exp.(-TE_vec/T2_vec[i])
    end
    
    return expSignal
    end
    
    T2_vec = 1:1:2000
    sDict = createExpBasis(Float32.(TEnum),Float32.(T2_vec))

    svd_obj = svd(sDict)
    basis = Complex.(svd_obj.V)[:,1:10]

    params = Dict{Symbol, Any}()
    params[:reconSize] = (N,N)
    params[:reco] = "multiCoilMultiEchoSubspace"
    #params[:sparseTrafo] = "nothing" #"nothing" #sparse trafo
    params[:reg] = L2Regularization(1.e-3)
    params[:iterations] = 1
    params[:solver] = CGNR
    params[:senseMaps] = reshape(coilsens,N,N,1,8)
    params[:basis] = basis

    α = reconstruction(acqData,params)
    x_approx2 = applySubspace(α,basis)

    relError = norm(x_approx - x_approx2)/norm(x_approx)
    @test relError < 5e-3
end

function testSpecificApplications(N=32)
    @testset "testSpecificApplications" begin
    testConversionAcqToRaw(N)
    testBinning(N)
    testSubspace(N)
    end
end

testSpecificApplications()
