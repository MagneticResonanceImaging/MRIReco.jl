@testset "Conversion" begin
  sx = 32
  sy = 24
  nSl = 6
  nEchos = 2
  T = ComplexF32
  kspace = ones(T,sx,sy,nSl,1,nEchos,1)

  # conversion 2D
  acq_2D = AcquisitionData(kspace,enc2D = true)
  @test acq_2D.traj[1].name == "Cartesian"
  @test size(acq_2D.kdata) == (nEchos,nSl,1)
  # conversion 3D
  acq_3D = AcquisitionData(kspace)
  @test acq_3D.traj[1].name == "Cartesian3D"
  @test size(acq_3D.kdata) == (nEchos,1,1)
end