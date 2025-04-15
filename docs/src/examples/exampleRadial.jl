using CairoMakie, MRIReco, MRISimulation
using ImageUtils: shepp_logan
include(joinpath(@__DIR__,"exampleUtils.jl"))
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

# use CairoMakie for  display
begin
  colormap=:grays
  f = Figure(size=(500,250))
  ax = Axis(f[1,1],title="Phantom")
  heatmap!(ax,rotr90(abs.(I));colormap)
  ax = Axis(f[1,2],title="Reconstruction")
  heatmap!(ax,rotr90(abs.(Ireco[:,:,1,1,1])))
  [hidedecorations!(f.content[ax]) for ax in eachindex(f.content)]
  f
end

# export images
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/simpleReco.png")
exportImage(filename, abs.(Ireco[:,:,1,1,1]))
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/phantom.png")
exportImage(filename, abs.(I[:,:,1,1,1]))
