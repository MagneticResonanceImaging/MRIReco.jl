
@testset "Flags" begin
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
  
  @test flag_is_set(p,"ACQ_IS_REVERSE") & flag_is_set(p,"ACQ_IS_NAVIGATION_DATA")
  flag_remove_all!(p)
  @test p.head.flags == UInt64(0)
  end
  