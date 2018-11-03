using MRIReco
using PyPlot

#### trajectories ####

trajectories = ["Spiral", "Radial", "Cartesian", "EPI", "SpiralVarDens"]
#                "Cartesian3D", "StackOfStars", "Kooshball"]

figure(1)
clf()

tr = trajectory("Spiral", 1, 600, windings=10)
subplot(2,2,1)
nodes = kspaceNodes(tr)
plot(nodes[1,:], nodes[2,:],"b.",lw=2)
title("Spiral")

tr = trajectory("Cartesian", 13, 50, EPI_factor=1)
subplot(2,2,2)
nodes = kspaceNodes(tr)
plot(nodes[1,:], nodes[2,:],"b.",lw=2)
title("Cartesian")


tr = trajectory("Radial", 13, 50)
subplot(2,2,3)
nodes = kspaceNodes(tr)
plot(nodes[1,:], nodes[2,:],"b.",lw=2)
title("Radial")


tr = trajectory("SpiralVarDens", 3, 200)
subplot(2,2,4)
nodes = kspaceNodes(tr)
plot(nodes[1,:], nodes[2,:],"b.",lw=2)
title("SpiralVarDens")



subplots_adjust(top=0.95,
       bottom=0.05,
       left=0.05,
       right=0.95,
       hspace=0.3,
       wspace=0.2)

savefig("../assets/trajectories.png", dpi=200)
