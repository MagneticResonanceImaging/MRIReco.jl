# Getting Started

We will start with a very simple example and perform simple simulation and
reconstruction based on a shepp logan phantom. The program looks like this
```julia
# image
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
```
We will go through the program step by step. First we create a 2D shepp logan
phantom of size `N=256`. Then we setup a dictionary that defines the simulation
parameters. Here, we chose a simple radial trajectory with 402 spokes and 512
samples per profile. We use a gridding-based simulator by setting `params[:simulation] = "fast"`

After setting up the parameter dictionary `params`, the simulation is performed
by calling
```julia
aqData = simulation(I, params)
```
The result `simulation` function outputs an acquisition object that is discussed
in more detail in the section [Acquisition Data](@ref).


### Reconstruction result
![Phantom](./assets/phantom.png)
![Reconstruction](./assets/simpleReco.png)
