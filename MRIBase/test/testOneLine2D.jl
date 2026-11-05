@testset "One-line density" begin
  for T in (Float32, Float64)
    for (numSamples, numProfiles, expectedProfile) in (
      (4, 3, [0.5, 0.25, 0.0, 0.25]),
      (5, 2, [0.5, 0.3, 0.1, 0.1, 0.3]),
      (1, 1, [0.5]),
    )
      density = oneLine2dDensity(T, numSamples, numProfiles)
      @test eltype(density) == T
      @test size(density) == (numSamples * numProfiles,)
      @test density ≈ repeat(T.(expectedProfile), numProfiles)
    end
  end
end
