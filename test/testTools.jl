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
  @testset "remove readout oversampling keeps the image in place" begin
    # a point simulated with 2x readout oversampling ends up where it is without
    # oversampling, at the center and at both edges of the readout direction
    N = 64
    params = Dict{Symbol, Any}(:reco => "direct", :reconSize => (N, N))
    for (x, y) in ((N ÷ 2 + 1, N ÷ 2 + 1), (1, 1), (N, 7))
      img = zeros(ComplexF64, N, N, 1)
      img[x, y, 1] = 1
      ref = reconstruction(simulation(trajectory(Float64, "Cartesian", N, N), img), params)
      imgOS = zeros(ComplexF64, 2N, N, 1)
      imgOS[x + N ÷ 2, y, 1] = 1
      raw = RawAcquisitionData(simulation(trajectory(Float64, "Cartesian", N, 2N), imgOS))
      raw.params["encodedFOV"] = [400.0, 200.0, 5.0]
      raw.params["reconFOV"] = [200.0, 200.0, 5.0]
      @test reconstruction(AcquisitionData(remove_oversampling(raw)), params) ≈ ref
    end
    # with equal fields of view, there is nothing to remove
    img = zeros(ComplexF64, N, N, 1)
    img[N ÷ 2 + 1, N ÷ 2 + 1, 1] = 1
    acq = simulation(trajectory(Float64, "Cartesian", N, N), img)
    raw = RawAcquisitionData(acq)
    raw.params["encodedFOV"] = [200.0, 200.0, 5.0]
    raw.params["reconFOV"] = [200.0, 200.0, 5.0]
    ref = reconstruction(acq, params)
    @test reconstruction(AcquisitionData(remove_oversampling(raw)), params) ≈ ref
  end
end