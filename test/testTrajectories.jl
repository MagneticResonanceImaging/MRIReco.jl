
@testset "Trajectories" begin

trajectories = ["Spiral", "Radial", "Cartesian", "EPI", "OneLine", "SpiralVarDens",
                "Cartesian3D", "StackOfStars", "Kooshball"]

for trajName in trajectories
  for T in [Float32, Float64]

    tr = trajectory(T, trajName, 8, 9)
    @test trajName == string(tr)
    #println(kspaceNodes(tr)) # TODO: make real tests
  end
end

end
