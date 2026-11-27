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
      pat(; kargs...) = sort(sample(shape, 3.0, pattern; kargs...))
      # by default, the default generator is used in its current state
      Random.seed!(3)
      idx = pat()
      @test pat() != idx
      Random.seed!(3)
      @test pat() == idx
      Random.seed!(3)
      @test pat(seed=nothing) == idx
    end
  end

  @testset "sample_kspace with $pattern" for pattern in ("random", "poisson", "vdPoisson"), profiles in (true, false)
    acqData = AcquisitionData(rand(ComplexF64, 64, 64, 1, 1, 2, 1))
    idx(; kargs...) = sample_kspace(acqData, 3.0, pattern; profiles=profiles, kargs...).subsampleIndices
    differ(pats) = pats[1] != pats[2]
    # a seed gives the same patterns, and each contrast gets its own pattern
    @test idx(seed=5) == idx(seed=5)
    @test idx(seed=5) != idx(seed=6)
    @test differ(idx(seed=5))
    @test !differ(idx(seed=5, rand=false))
    # the same with a generator, which advances
    @test idx(rng=StableRNG(1)) == idx(rng=StableRNG(1))
    @test idx(rng=StableRNG(1)) != idx(rng=StableRNG(2))
    @test differ(idx(rng=StableRNG(1)))
    @test !differ(idx(rng=StableRNG(1), rand=false))
    rng = StableRNG(1)
    @test idx(rng=rng, rand=false) != idx(rng=rng, rand=false)
    @test idx(rng=StableRNG(1), seed=5) == idx(rng=StableRNG(2), seed=5)
    # by default, the default generator is used in its current state
    Random.seed!(3)
    pats = idx()
    @test differ(pats)
    @test idx() != pats
    @test !differ(idx(rand=false))
    Random.seed!(3)
    @test idx() == pats
    # the same in place
    inPlace = deepcopy(acqData)
    Random.seed!(3)
    MRISampling.sample_kspace!(inPlace, 3.0, pattern; profiles=profiles)
    @test inPlace.subsampleIndices == pats
  end

  @testset "Variable density Poisson disk $pattern" for (pattern, opts) in
      (("vdPoissonBART", (vardens=5.0,)), ("vdPoissonGrid", (vardens=2.0,)))
    for shape in ((256,), (64,64))
      pat(; kargs...) = sort(sample(shape, 3.0, pattern; opts..., kargs...))
      # seed and generator, as for the other random patterns
      idx = pat(seed=1234)
      @test pat(seed=1234) == idx
      @test pat(seed=4321) != idx
      @test pat(rng=StableRNG(1)) == pat(rng=StableRNG(1))
      @test pat(rng=StableRNG(1)) != pat(rng=StableRNG(2))
      rng = StableRNG(1)
      @test pat(rng=rng) != pat(rng=rng)
      @test pat(rng=StableRNG(1), seed=5) == pat(rng=StableRNG(2), seed=5)
      Random.seed!(1)
      expected = rand()
      Random.seed!(1)
      pat(rng=StableRNG(1))
      @test rand() == expected
      # exact number of samples, no duplicates
      @test length(idx) == floor(Int, prod(shape) / 3.0)
      @test allunique(idx)
    end

    # calibration area, ellipse and variable density
    M = 64
    idx = sample((M,M), 4.0, pattern; calsize=8, elliptical=true, seed=1, opts...)
    mask = falses(M,M)
    mask[idx] .= true
    @test length(idx) == floor(Int, M^2 / 4.0)
    @test all(mask[29:36, 29:36])
    ρ = [sqrt(((y-33)/32)^2 + ((z-33)/32)^2) for y=1:M, z=1:M]
    # vdPoissonBART crops the ellipse before rounding down to the grid, as BART: the grid point y covers the offsets
    # [y-33, y-32) from the center, and its offset closest to the center must be inside the ellipse
    near(y) = y >= 33 ? y - 33 : 32 - y
    ρnear = [sqrt((near(y)/32)^2 + (near(z)/32)^2) for y=1:M, z=1:M]
    @test !any(mask[(pattern == "vdPoissonBART" ? ρnear : ρ) .> 1])
    density(sel) = count(mask[sel]) / count(sel)
    @test density(0.15 .< ρ .< 0.4) > density(0.6 .< ρ .<= 1)

    # the center of k-space is always sampled, also without calibration area
    @test all(129 in sample((256,), 3.0, pattern; seed=s, opts...) for s in 1:20)
    @test all(LinearIndices((M,M))[33,33] in sample((M,M), 4.0, pattern; seed=s, opts...) for s in 1:20)

    # too small reduction factor, 3D shape
    strong = (vardens = pattern == "vdPoissonBART" ? 40.0 : 6.0,)
    @test_throws ArgumentError sample((M,M), 1.5, pattern; strong...)
    @test_throws ErrorException sample((M,M,M), 4.0, pattern)

    # 3D acquisition: the pattern of the phase encoding plane, for all readout samples
    idx3 = sample_lines((16,M,M), 4.0; sampleFunc=pattern, calsize=8, seed=1, opts...)
    @test length(idx3) == 16 * floor(Int, M^2 / 4.0)

    # undersampling of acquisition data, with profiles or single samples
    acqData = AcquisitionData(rand(ComplexF64, 64, 64, 1, 1, 2, 1))
    for profiles in (true, false)
      sub(; kargs...) = sample_kspace(acqData, 3.0, pattern; profiles=profiles, opts..., kargs...).subsampleIndices
      pats = sub(seed=5)
      @test sub(seed=5) == pats
      @test pats[1] != pats[2]
      @test sub(seed=5, rand=false)[1] == sub(seed=5, rand=false)[2]
    end
  end

  @testset "Patterns without random numbers" begin
    acqData = AcquisitionData(rand(ComplexF64, 64, 64, 1, 1, 2, 1))
    idx = sample_kspace(acqData, 2.0, "regular").subsampleIndices
    @test sample_kspace(acqData, 2.0, "regular"; rng=StableRNG(1)).subsampleIndices == idx
    @test sample_kspace(acqData, 2.0, "regular"; seed=5).subsampleIndices == idx
  end
end
