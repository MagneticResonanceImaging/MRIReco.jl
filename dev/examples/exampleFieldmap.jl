using MRIReco

##### simple example ####

N = 256
I = shepp_logan(N)
I = circularShutterFreq!(I,1)
cmap = 1im*quadraticFieldmap(N,N,125*2pi)

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Spiral"
params[:numProfiles] = 1
params[:numSamplingPerProfile] = N*N
params[:windings] = 128
params[:AQ] = 3.0e-2
params[:correctionMap] = cmap[:,:,1]

# do simulation
acqData = simulation(I, params)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N,N)
Ireco = reconstruction(acqData, params)

# export image
Icolored = colorview(Gray, abs.(Ireco)./maximum(abs.(Ireco)))
save("../assets/fieldmapReco1.png", Icolored[:,:,1,1,1] )

# now with fieldmap correction
params[:reco] = "standard"
params[:regularization] = "L2"
params[:iterations] = 3
params[:solver] = "admm"
params[:correctionMap] = cmap
params[:alpha] = 1.75
params[:m] = 4.0
params[:K] = 28

Ireco = reconstruction(acqData, params)


# export images
Icolored = colorview(Gray, abs.(Ireco)./maximum(abs.(Ireco)))
save("../assets/fieldmapReco2.png", Icolored[:,:,1,1,1] )


# export fieldmap image
cmap_ = imag(cmap[:,:,1]) .+ 125*2pi
cmapColored = colorview(Gray, abs.(cmap_)./maximum(abs.(cmap_)))
save("../assets/fieldmap.png", cmapColored )
