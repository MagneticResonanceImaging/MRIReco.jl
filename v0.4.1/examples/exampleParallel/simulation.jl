using PyPlot, MRIReco, FFTW, LinearAlgebra, BenchmarkTools

FFTW.set_num_threads(1);BLAS.set_num_threads(1)

N = 256
numCoils = 8
I = shepp_logan(N)
I = circularShutterFreq!(I,1)

coilsens = birdcageSensitivity(N, numCoils, 1.5)
I = circularShutterFreq!(I,1)
cmap = 1im*quadraticFieldmap(N,N,125*2pi)

# simulation parameters
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Spiral"
params[:numProfiles] = 6
params[:numSamplingPerProfile] = div(N*N,16)
params[:windings] = div(N,16)
params[:AQ] = 2.0e-2
params[:senseMaps] = coilsens
params[:verbose] = false
#params[:correctionMap] = cmap #[:,:,1]

# do simulation
acqData = simulation(I, params)

raw = RawAcquisitionData(acqData)
fout = ISMRMRDFile(joinpath(@__DIR__(),"simparallel.h5"))
save(fout, raw)

