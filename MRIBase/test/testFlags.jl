
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
  end

