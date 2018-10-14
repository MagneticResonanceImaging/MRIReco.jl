# Getting Started

We will start with a very simple example and perform simple simulation and
reconstruction based on a shepp logan phantom. The programm looks like this
```julia
# image
N = 256
x = shepp_logan(N)

# simulation
params = Dict{Symbol, Any}()
params[:simulation] = "nfft"
params[:trajName] = "Cartesian"
params[:numProfiles] = floor(Int64, N)
params[:numSamplingPerProfile] = N

aqData = simulation(x, params)
aqData = MRIReco.sample_kspace(aqData, redFac, "poisson", calsize=5)
aqData = convertUndersampledData(aqData)

# reco
params[:reco] = "simple"    # encoding model
params[:shape] = (N,N)
params[:sparseTrafoName] = "Wavelet" #sparse trafo
params[:regularization] = "L1"       # regularization
params[:lambdL1] = 1.e-3
params[:solver] = "admm"    # solver
params[:iterations] = 1000
params[:ρ] = 1.0e-1

x_reco = reconstruction(aqData, params)
```
We will go through the program step by step. First we create a 2D shepp logan
phantom of size N=256. Then we setup a dictionary that defines the simulation
parameters. Here, we chose a simple Cartesian trajectory ...
