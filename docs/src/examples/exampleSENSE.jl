using MRIReco

##### simple example ####

N = 256
numCoils = 8
I = shepp_logan(N)
I = circularShutterFreq!(I,1)
cmap = 1im*quadraticFieldmap(N,N,125*2pi)

coilsens = birdcageSensitivity(N, 8, 1.5)

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Spiral"
params[:numProfiles] = 1
params[:numSamplingPerProfile] = div(N*N,2)
params[:windings] = div(N,4)
params[:AQ] = 3.0e-2
#params[:correctionMap] = cmap[:,:,1]
params[:senseMaps] = coilsens

# do simulation
acqData = simulation(I, params)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil" #"standard"
params[:shape] = (N,N)
params[:regularization] = "L2"
params[:iterations] = 10
params[:solver] = "admm"
#params[:correctionMap] = cmap
params[:senseMaps] = reshape(coilsens, N*N, numCoils, 1)
params[:alpha] = 1.75
params[:m] = 4.0
params[:K] = 28
Ireco = reconstruction(acqData, params)

# export images
Icolored = colorview(Gray, abs.(Ireco)./maximum(abs.(Ireco)))
save("../assets/senseReco.png", Icolored[:,:,1,1,1] )
