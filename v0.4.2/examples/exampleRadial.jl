using PyPlot, MRIReco 

##### simple example ####

N = 256
I = shepp_logan(N)

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Radial"
params[:numProfiles] = floor(Int64, pi/2*N)
params[:numSamplingPerProfile] = 2*N

# do simulation
acqData = simulation(I, params)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N,N)
Ireco = reconstruction(acqData, params)

# use PyPlot for interactive display
figure(1)
clf()
subplot(1,2,1)
title("Phantom")
imshow(abs.(I))
subplot(1,2,2)
title("Reconstruction")
imshow(abs.(Ireco[:,:,1,1,1]))

# export images
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/simpleReco.png")
exportImage(filename, abs.(Ireco[:,:,1,1,1]))
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/phantom.png")
exportImage(filename, abs.(I[:,:,1,1,1]))
