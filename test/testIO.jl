
@testset "IO" begin

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
params[:alpha] = 1.75
params[:m] = 4.0
params[:K] = 28
Ireco = reconstruction(acqData, params)

saveImage("testReco.nii", Ireco, false)
IrecoLoaded = loadImage("testReco.nii")
@test (norm(vec(Ireco)-vec(IrecoLoaded))/norm(vec(Ireco))) < 1e-5

saveImage("testRecoAbs.nii", Ireco, true)
IrecoLoadedAbs = loadImage("testRecoAbs.nii")
@test (norm(vec(abs.(Ireco))-vec(IrecoLoadedAbs))/norm(vec(abs.(Ireco)))) < 1e-5

end
