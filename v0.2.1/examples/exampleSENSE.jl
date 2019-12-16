using PyPlot, MRIReco

##### simple example ####

N = 256
numCoils = 8
I = shepp_logan(N)
I = circularShutterFreq!(I,1)

coilsens = birdcageSensitivity(N, 8, 1.5)

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Spiral"
params[:numProfiles] = 6
params[:numSamplingPerProfile] = div(N*N,16)
params[:windings] = div(N,16)
params[:AQ] = 2.0e-2
params[:senseMaps] = coilsens

# do simulation
acqData = simulation(I, params)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:λ] = 1.e-3
params[:iterations] = 40
params[:solver] = "cgnr"
params[:senseMaps] = coilsens
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
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/senseReco.png")
exportImage(filename, abs.(Ireco[:,:,1,1,1]) )
