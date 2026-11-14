using Test, MRIBase, MRICoilSensitivities, MRISampling, MRISimulation, MRIReco
using LinearAlgebra, ImageUtils
using MRICoilSensitivities.FFTW: fft, ifft, fftshift, ifftshift

include("testCoilCompression.jl")
include("testCoilsens.jl")
