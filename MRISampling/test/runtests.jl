using MRISampling, MRIBase, Test

@testset "MRISampling" begin
  @testset "Seed of $pattern" for pattern in ("random", "poisson", "vdPoisson")
    for shape in ((256,), (64,64))
      idx = sort(sample(shape, 3.0, pattern; seed=1234))
      @test sort(sample(shape, 3.0, pattern; seed=1234)) == idx
      @test sort(sample(shape, 3.0, pattern; seed=4321)) != idx
    end
  end

  @testset "Seed of sample_kspace with $pattern" for pattern in ("random", "poisson", "vdPoisson")
    acqData = AcquisitionData(rand(ComplexF64, 64, 64, 1, 1, 2, 1))
    # the default seed gives the same patterns, and each contrast gets its own pattern
    idx = sample_kspace(acqData, 3.0, pattern).subsampleIndices
    @test sample_kspace(acqData, 3.0, pattern).subsampleIndices == idx
    @test idx[1] != idx[2]
  end
end
