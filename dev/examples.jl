using MRIReco
using ImageView

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
aqData = simulation(I, params)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "nfft"
params[:shape] = (N,N)
params[:alpha] = 1.75
params[:m] = 4.0
params[:K] = 28
Ireco = reconstruction(aqData, params)

# use ImageView for interactive display
imshow(abs.(I))
imshow(abs.(Ireco[:,:,1,1,1]))

# export images
Icolored = colorview(Gray, abs.(Ireco)./maximum(abs.(Ireco)))
save("assets/simpleReco.png", Icolored[:,:,1,1,1] )

Icolored = colorview(Gray, abs.(I)./maximum(abs.(I)))
save("assets/phantom.png", Icolored[:,:] )

#### trajectories ####
