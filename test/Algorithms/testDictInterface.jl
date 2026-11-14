# The dict interface `reconstruction(acqData, ::Dict)` translates into reco plans. These tests cover
# translations that the plan tests do not reach.

function testDictSparseTrafoOperator(N = 32)
  x = ComplexF64.(shepp_logan(N))
  params = Dict{Symbol, Any}(:simulation => "fast", :trajName => "Cartesian", :numProfiles => N, :numSamplingPerProfile => N)
  acqData = simulation(real(x), params)

  params = Dict{Symbol, Any}(:reco => "standard", :reg => L1Regularization(5.e-2), :solver => FISTA, :iterations => 20)
  Inone = reconstruction(acqData, params)
  Istring = reconstruction(acqData, merge(params, Dict{Symbol, Any}(:sparseTrafo => "Wavelet")))
  Iop = reconstruction(acqData, merge(params, Dict{Symbol, Any}(:sparseTrafo => WaveletOp(ComplexF64; shape = (N, N)))))

  @test Iop ≈ Istring
  @test norm(vec(Iop) - vec(Inone)) / norm(vec(Inone)) > 1e-2
end

function testDictConverted3d(N = 32, Nz = 8)
  x = ComplexF64.(shepp_logan(N))
  I = permutedims(cat([(1 - 0.1 * i) .* x for i in 0:Nz-1]..., dims = 3), [3, 1, 2])
  senseMaps = zeros(ComplexF64, Nz, N, N, 2)
  senseMaps[:, 1:N÷2, :, 1] .= 1
  senseMaps[:, N÷2+1:end, :, 2] .= 1

  params = Dict{Symbol, Any}(:simulation => "fast", :trajName => "Cartesian3D", :numProfiles => N,
                             :numSamplingPerProfile => Nz, :numSlices => N, :senseMaps => senseMaps)
  acqData = convert3dTo2d(simulation(real(I), params))

  # reconSize and sensitivity maps given in 3D order
  params = Dict{Symbol, Any}(:reco => "multiCoil", :reconSize => (Nz, N, N), :senseMaps => senseMaps,
                             :reg => L2Regularization(1.e-3), :solver => CGNR, :iterations => 10)
  Ireco = reconstruction(acqData, params)

  @test size(Ireco)[1:3] == (N, N, Nz)
  Iref = permutedims(I, [2, 3, 1])
  @test norm(vec(abs.(Ireco)) - vec(abs.(Iref))) / norm(vec(Iref)) < 1e-2
end

function testReconstructByName(N = 32)
  x = ComplexF64.(shepp_logan(N))
  params = Dict{Symbol, Any}(:simulation => "fast", :trajName => "Cartesian", :numProfiles => N, :numSamplingPerProfile => N)
  acqData = simulation(real(x), params)

  reg = L2Regularization(1.e-3)
  Idict = reconstruction(acqData, Dict{Symbol, Any}(:reco => "standard", :reg => reg, :solver => CGNR, :iterations => 5))
  @test reconstruct("standard", acqData; reg, solver = CGNR, iterations = 5) ≈ Idict

  planfile = joinpath(pkgdir(MRIReco), "reco_plans", "standard.toml")
  @test MRIRecoPlan(planfile) isa RecoPlan{StandardReconstruction}
  @test MRIRecoPlan("standard") isa RecoPlan{StandardReconstruction}
end

function testDictNoWarnings(N = 32)
  x = ComplexF64.(shepp_logan(N))
  params = Dict{Symbol, Any}(:simulation => "fast", :trajName => "Cartesian", :numProfiles => N, :numSamplingPerProfile => N)
  acqData = simulation(real(x), params)

  # defaultRecoParams sets solver parameters, which a direct reconstruction does not use
  @test_logs min_level = Base.CoreLogging.Warn reconstruction(acqData, Dict{Symbol, Any}(:reco => "direct"))
  @test_logs min_level = Base.CoreLogging.Warn reconstruction(acqData, Dict{Symbol, Any}(:reco => "standard", :iterations => 2, :unknownKey => 1))
end

@testset "Dict interface" begin
  @testset "sparseTrafo as operator" testDictSparseTrafoOperator()
  @testset "3D data converted to 2D" testDictConverted3d()
  @testset "reconstruct by plan name" testReconstructByName()
  @testset "no warnings for unused keys" testDictNoWarnings()
end
