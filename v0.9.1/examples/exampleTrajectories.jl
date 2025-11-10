using CairoMakie
using MRIReco

#### trajectories ####

trajectories = ["Spiral", "Radial", "Cartesian", "EPI", "SpiralVarDens"]

f = Figure()

T=Float32

tr = trajectory(Float32,"Spiral", 1, 600, windings=10)
nodes = kspaceNodes(tr)
ax=Axis(f[1,1],title="Spiral")
plot!(ax,nodes[1,:], nodes[2,:])


tr = trajectory(T,"Cartesian", 13, 50, EPI_factor=1)
nodes = kspaceNodes(tr)
ax=Axis(f[1,2],title="Cartesian")
plot!(ax,nodes[1,:], nodes[2,:])



tr = trajectory(T,"Radial", 13, 50)
nodes = kspaceNodes(tr)
ax=Axis(f[2,1],title="Radial")
plot!(ax,nodes[1,:], nodes[2,:])




tr = trajectory(T,"SpiralVarDens", 3, 200)
nodes = kspaceNodes(tr)
ax=Axis(f[2,2],title="SpiralVarDens")
plot!(ax,nodes[1,:], nodes[2,:])
f


filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/trajectories.png")
save(filename,f)
