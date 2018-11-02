# Parallel Imaging

For parallel imaging MRIReco.jl uses an iterative SENSE approach. In the following
code example we show how to simulate MRI data with an array of 8 coils and a corresponding
SENSE reconstruction.

```julia
N = 256
numCoils = 8
I = shepp_logan(N)
I = circularShutterFreq!(I,1)

coilsens = birdcageSensitivity(N, 8, 1.5)

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Spiral"
params[:numProfiles] = 1
params[:numSamplingPerProfile] = div(N*N,2)
params[:windings] = div(N,4)
params[:AQ] = 3.0e-2
params[:senseMaps] = coilsens

# do simulation
acqData = simulation(I, params)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:shape] = (N,N)
params[:regularization] = "L2"
params[:iterations] = 10
params[:solver] = "admm"
params[:senseMaps] = coilsens

Ireco = reconstruction(acqData, params)

```
