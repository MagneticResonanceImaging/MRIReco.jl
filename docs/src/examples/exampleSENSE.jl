using CairoMakie, MRIReco, MRIReco.RegularizedLeastSquares
using MRICoilSensitivities, MRISimulation

include(joinpath(@__DIR__,"exampleUtils.jl"))
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
params[:reg] = L2Regularization(1.e-1)
params[:iterations] = 1
params[:solver] = CGNR
params[:senseMaps] = coilsens
params[:toeplitz] = false
params[:normalizeReg] = MeasurementBasedNormalization()

@time Ireco_u = reconstruction(acqData, params)

params[:iterations] = 40
@time Ireco = reconstruction(acqData, params)

# use CairoMakie for  display
begin
  colormap=:grays
  f = Figure(size=(750,250))
  ax = Axis(f[1,1],title="Phantom")
  heatmap!(ax,rotr90(abs.(I));colormap)
  ax = Axis(f[1,2],title="Reconstruction iteration 1")
  heatmap!(ax,rotr90(abs.(Ireco_u[:,:,1,1,1]));colormap)
  ax = Axis(f[1,3],title="Reconstruction iteration 40")
  heatmap!(ax,rotr90(abs.(Ireco[:,:,1,1,1]));colormap)
  [hidedecorations!(f.content[ax]) for ax in eachindex(f.content)]
  f
end

# export images
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/senseReco.png")
exportImage(filename, abs.(Ireco[:,:,1,1,1]) )
