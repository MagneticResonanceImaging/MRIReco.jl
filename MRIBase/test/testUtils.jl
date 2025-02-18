
@testset "filter_raw_by_flags" begin

  head = AcquisitionHeader()
  traj = Matrix{Float32}(undef,0,0)
  dat =  Matrix{ComplexF32}(undef,0,0)

  # Test data
  flags = [
    ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING | ACQ_IS_REVERSE,
    ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING,
    ACQ_IS_REVERSE
  ]

  profiles = [Profile(AcquisitionHeader(flags = flag), traj, dat) for flag in flags]

  rawData = RawAcquisitionData(minimalHeader((128,128),(250.,250.,250.)), profiles)

  # Test case 1: Filter by one flag
  filteredData = filter_raw_by_flags(rawData, ["ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING"])
  @test length(filteredData.profiles) == 2
 
  @test ~all(profile -> flag_is_set(profile,"ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING"),rawData.profiles)
  @test all(profile -> flag_is_set(profile,"ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING"),filteredData.profiles)

  # Test case 2: Filter by two flags
  filteredData = filter_raw_by_flags(rawData, ["ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING", "ACQ_IS_REVERSE"])
  @test length(filteredData.profiles) == 1
  @test all(profile -> (flag_is_set(profile,"ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING") &&flag_is_set(profile,"ACQ_IS_REVERSE")),filteredData.profiles)

  # Test case 3: Filter by a flag that does not exist
  filteredData = filter_raw_by_flags(rawData, ["ACQ_USER8"])
  @test length(filteredData.profiles) == 0
  
  # Test case 4: Filter by a flag name  that does not exist
  @test_throws Exception filteredData = filter_raw_by_flags(rawData, "NON_EXISTENT_FLAG")
end