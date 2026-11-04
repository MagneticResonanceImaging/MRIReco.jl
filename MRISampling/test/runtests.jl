using MRISampling, MRIBase, Random, StableRNGs, Test

@testset "MRISampling" begin
  @testset "Seed of $pattern" for pattern in ("random", "poisson", "vdPoisson")
    for shape in ((256,), (64,64))
      idx = sort(sample(shape, 3.0, pattern; seed=1234))
      @test sort(sample(shape, 3.0, pattern; seed=1234)) == idx
      @test sort(sample(shape, 3.0, pattern; seed=4321)) != idx
    end
  end

  @testset "Generator of $pattern" for pattern in ("random", "poisson", "vdPoisson")
    for shape in ((256,), (64,64))
      pat(; kargs...) = sort(sample(shape, 3.0, pattern; kargs...))
      idx = pat(rng=StableRNG(1))
      @test pat(rng=StableRNG(1)) == idx
      # a generator isn't seeded by default, and its state advances
      @test pat(rng=StableRNG(2)) != idx
      rng = StableRNG(1)
      @test pat(rng=rng) == idx
      @test pat(rng=rng) != idx
      # a seed sets its state
      @test pat(rng=StableRNG(1), seed=5) == pat(rng=StableRNG(2), seed=5)
      @test pat(rng=StableRNG(1), seed=5) != pat(rng=StableRNG(1), seed=6)
      # the default generator isn't touched
      Random.seed!(1)
      expected = rand()
      Random.seed!(1)
      pat(rng=StableRNG(1))
      pat(rng=StableRNG(1), seed=5)
      @test rand() == expected
    end
  end

  @testset "Generator of the other methods" begin
    lines(rng) = sample_lines((16,32,8), 3.0; sampleFunc="poisson", rng=rng)
    @test lines(StableRNG(1)) == lines(StableRNG(1))
    @test lines(StableRNG(1)) != lines(StableRNG(2))
    array(rng) = sample_kspace(zeros(ComplexF64, 64, 64), 3.0, "poisson"; rng=rng)[2]
    @test array(StableRNG(1)) == array(StableRNG(1))
    @test array(StableRNG(1)) != array(StableRNG(2))
  end

  @testset "No seed for $pattern" for pattern in ("random", "poisson", "vdPoisson")
    for shape in ((256,), (64,64))
      pat() = sort(sample(shape, 3.0, pattern; seed=nothing))
      Random.seed!(3)
      idx = pat()
      @test pat() != idx
      Random.seed!(3)
      @test pat() == idx
    end
  end

  @testset "Default seed" begin
    @test sample((64,64), 3.0, "poisson") == sample((64,64), 3.0, "poisson")
    @test sample((64,64), 3.0, "vdPoisson") == sample((64,64), 3.0, "vdPoisson")
    @test sample((64,64), 3.0, "random") != sample((64,64), 3.0, "random")
  end

  @testset "sample_kspace with $pattern" for pattern in ("random", "poisson", "vdPoisson"), profiles in (true, false)
    acqData = AcquisitionData(rand(ComplexF64, 64, 64, 1, 1, 2, 1))
    idx(; kargs...) = sample_kspace(acqData, 3.0, pattern; profiles=profiles, kargs...).subsampleIndices
    differ(pats) = pats[1] != pats[2]
    # the default seed gives the same patterns, and each contrast gets its own pattern
    @test idx() == idx()
    @test differ(idx())
    @test idx(seed=5) == idx(seed=5)
    @test idx(seed=5) != idx(seed=6)
    @test !differ(idx(seed=5, rand=false))
    # the same with a generator, which isn't seeded by default and advances
    @test idx(rng=StableRNG(1)) == idx(rng=StableRNG(1))
    @test idx(rng=StableRNG(1)) != idx(rng=StableRNG(2))
    @test differ(idx(rng=StableRNG(1)))
    @test !differ(idx(rng=StableRNG(1), rand=false))
    rng = StableRNG(1)
    @test idx(rng=rng, rand=false) != idx(rng=rng, rand=false)
    @test idx(rng=StableRNG(1), seed=5) == idx(rng=StableRNG(2), seed=5)
    # the same without a seed for the default generator
    @test idx(seed=nothing) != idx(seed=nothing)
    @test differ(idx(seed=nothing))
    @test !differ(idx(seed=nothing, rand=false))
    Random.seed!(3)
    first = idx(seed=nothing)
    Random.seed!(3)
    @test idx(seed=nothing) == first
  end

  @testset "Patterns without random numbers" begin
    acqData = AcquisitionData(rand(ComplexF64, 64, 64, 1, 1, 2, 1))
    idx = sample_kspace(acqData, 2.0, "regular").subsampleIndices
    @test sample_kspace(acqData, 2.0, "regular"; rng=StableRNG(1)).subsampleIndices == idx
    @test sample_kspace(acqData, 2.0, "regular"; seed=nothing).subsampleIndices == idx
  end
end
