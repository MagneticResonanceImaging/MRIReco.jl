
@testset "Trajectories" begin

trajectories = ["Spiral", "Radial", "Cartesian", "EPI", "OneLine", "SpiralVarDens",
                "Cartesian3D", "StackOfStars", "Kooshball"]

for trajName in trajectories

  tr = trajectory(trajName, 8, 9)
  @test trajName == string(tr)
  #println(kspaceNodes(tr)) # TODO: make real tests

end

end
