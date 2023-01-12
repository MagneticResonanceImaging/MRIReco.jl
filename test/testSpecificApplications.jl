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

function testSpecificApplications(N=32)
    @testset "testSpecificApplications" begin
    testConversionAcqToRaw(N)
    testBinning(N)
    end
end

testSpecificApplications()
