using MRIBase
@testset "Tools" begin
  @testset "remove readout oversampling" begin
    NRead = 128
    Oversampling = 2
    nChan = 4
    head = AcquisitionHeader()
    traj = Matrix{Float32}(undef,0,0)
    dat =  ones(ComplexF32,NRead*Oversampling,nChan)
  
    p = Profile(head,traj,dat)
    params = minimalHeader((NRead*Oversampling,96),(500.,250.,250.))
    params["reconFOV"] = (250.,250.,250.)
    rawData = RawAcquisitionData(params, [p])
    rawData2 = remove_oversampling(rawData)
    @test size(rawData2.profiles[1].data,1) == NRead
  end
end