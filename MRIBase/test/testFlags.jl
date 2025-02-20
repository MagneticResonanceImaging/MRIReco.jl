
@testset "Flags" begin
  @testset "Flags one by one" begin
  head = AcquisitionHeader()
  traj = Matrix{Float32}(undef,0,0)
  dat =  Matrix{ComplexF32}(undef,0,0)
  p = Profile(head,traj,dat)

  p2 = deepcopy(p)
  @test ~flag_is_set(p,"ACQ_IS_REVERSE")
  flag_set!(p,"ACQ_IS_REVERSE")
  flag_set!(p2,22)
  @test flag_is_set(p,"ACQ_IS_REVERSE")
  @test p2.head.flags == p.head.flags
  flag_remove!(p,"ACQ_IS_REVERSE")
  @test ~flag_is_set(p,22)

  flag_set!(p,"ACQ_IS_REVERSE")
  flag_set!(p,"ACQ_IS_NAVIGATION_DATA")

  fl = flags_of(p)

  @test isequal(fl[1],"ACQ_IS_NAVIGATION_DATA") && isequal(fl[2],"ACQ_IS_REVERSE")
  @test flag_is_set(p,"ACQ_IS_REVERSE") & flag_is_set(p,"ACQ_IS_NAVIGATION_DATA")
  flag_remove_all!(p)
  @test p.head.flags == UInt64(0)
  end

  @testset "Flags vector" begin
    head = AcquisitionHeader()
    traj = Matrix{Float32}(undef,0,0)
    dat =  Matrix{ComplexF32}(undef,0,0)
    p = Profile(head,traj,dat)
    p2 = deepcopy(p)
    p3 = deepcopy(p2)

    flag_set!(p,["ACQ_IS_REVERSE","ACQ_USER8"])

    flag_set!(p2,"ACQ_IS_REVERSE")
    flag_set!(p2,"ACQ_USER8")

    flag_set!(p3,[FLAGS["ACQ_USER8"],FLAGS["ACQ_IS_REVERSE"]])
    @test p.head.flags == p2.head.flags == p3.head.flags

    flag_remove!(p,["ACQ_IS_REVERSE"])
    flag_remove!(p2,"ACQ_IS_REVERSE")
    flag_remove!(p3,FLAGS["ACQ_IS_REVERSE"])
    @test p.head.flags == p2.head.flags == p3.head.flags

     flags_of(p)
  flag_remove!(p,FLAGS["ACQ_USER8"])
  flags_of(p)
  end

  @testset "Flags Const interface" begin
    head = AcquisitionHeader()
    traj = Matrix{Float32}(undef,0,0)
    dat =  Matrix{ComplexF32}(undef,0,0)
    p = Profile(head,traj,dat)
    p2 = deepcopy(p)

    flag_set!(p,["ACQ_IS_REVERSE","ACQ_USER8"])
    p2.head.flags += ACQ_IS_REVERSE + ACQ_USER8
    @test p.head.flags == p2.head.flags

    flag_remove!(p,"ACQ_IS_REVERSE")
    p2.head.flags += -ACQ_IS_REVERSE
    @test p.head.flags == p2.head.flags
  end

  @testset "Test all flags" begin
    for (key, val) in FLAGS
      # Initialize profiles
      head = AcquisitionHeader()
      traj = Matrix{Float32}(undef,0,0)
      dat =  Matrix{ComplexF32}(undef,0,0)
      p1 = Profile(head,traj,dat)
      p2 = deepcopy(p1) 
      # Setting flags
      flag_set!(p1, val)
      p2.head.flags = 1b64 << ( val - 1 )
      @test p1.head.flags == p2.head.flags
  end
  end

  @testset "Flags errors" begin
    head = AcquisitionHeader()
    traj = Matrix{Float32}(undef,0,0)
    dat =  Matrix{ComplexF32}(undef,0,0)
    p = Profile(head,traj,dat)
    @test_throws Exception flag_set!(p, -1)
    @test_throws Exception flag_set!(p, "NOT_A_VALID_FLAG")
    @test_throws Exception flag_is_set(p, -1)
    @test_throws Exception flag_is_set(p, "NOT_A_VALID_FLAG")
    @test_throws Exception flag_remove!(p, -1)
    @test_throws Exception flag_remove!(p, "NOT_A_VALID_FLAG")
    @test_throws Exception flag_remove!(p,UInt64(0))
    @test_throws Exception flag_remove!(p,ACQ_IS_REVERSE)
    @test_throws DomainError p.head.flags = "test"
  end

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
    
    @testset "remove_raw_by_flags" begin
      head = AcquisitionHeader()
      traj = Matrix{Float32}(undef,0,0)
      dat =  Matrix{ComplexF32}(undef,0,0)
    
      # Test data
      flags = [
        ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING | ACQ_IS_REVERSE,
        ACQ_IS_NOISE_MEASUREMENT,
        ACQ_IS_REVERSE
      ]
      profiles = [Profile(AcquisitionHeader(flags = flag), traj, dat) for flag in flags]
      rawData = RawAcquisitionData(minimalHeader((128,128),(250.,250.,250.)), profiles)
    
      filteredData = remove_raw_by_flags(rawData)
      @test length(filteredData.profiles) == 2
      filteredData = remove_raw_by_flags(rawData,[])
      @test length(filteredData.profiles) == 3
      filteredData = remove_raw_by_flags(rawData,["ACQ_IS_REVERSE"])
      @test length(filteredData.profiles) == 1
      filteredData = remove_raw_by_flags(rawData,"ACQ_IS_REVERSE")
      @test length(filteredData.profiles) == 1
    end
  end

