using MRIReco, FFTW, LinearAlgebra, BenchmarkTools

FFTW.set_num_threads(1);BLAS.set_num_threads(1)


N = 256
numCoils = 8
coilsens = birdcageSensitivity(N, numCoils, 1.5)
cmap = 1im*quadraticFieldmap(N,N,125*2pi)

f = ISMRMRDFile(joinpath(@__DIR__(),"simparallel.h5"))
raw = RawAcquisitionData(f)
acqData = AcquisitionData(raw)

# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:λ] = 1.e-3
params[:iterations] = 40
params[:solver] = "cgnr"
params[:senseMaps] = coilsens
#params[:correctionMap] = cmap

t = @belapsed global img = reconstruction(acqData, params)
filename = joinpath(@__DIR__(),"reco.png")
exportImage(filename, abs.(img[:,:,1,1,1]) )

println(t)


