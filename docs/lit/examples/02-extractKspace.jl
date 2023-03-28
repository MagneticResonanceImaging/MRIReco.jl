#---------------------------------------------------------
# # [Extract k-space](@id 01-extractKspace)
#---------------------------------------------------------


#=
## Description

This example described how to extract kspace from cartesian datasets,
undersample the data and create back and AcquisitionData structure that can be used for reconstruction.

This possibility can be easily combined with the [BartIO](https://github.com/MagneticResonanceImaging/BartIO.jl) package to reconstruct the kspace.
=#

# ## Setup
using CairoMakie
using ImageUtils: shepp_logan
using MRIReco, MRISimulation
using InteractiveUtils: versioninfo

function plot_im2D(im2D;title::String="")
    f = Figure()
    ax = Axis(f[1, 1],aspect = DataAspect(), yreversed = true, title = title)
    image!(ax, im2D')
    hidedecorations!(ax, grid = false)
    f
end

# Let's create a non-square AcquisitionData structure and perform a standard reconstruction
N = 128
N2 = 96

# image
x = shepp_logan(N)

# simulation
params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Cartesian"
params[:numSamplingPerProfile] = N
params[:numProfiles] = floor(Int64, N2)

acqData = simulation(x[1:N,1:N2], params)

params[:reco] = "direct"
params[:reconSize] = (N,N2)

x_approx = reconstruction(acqData, params);

# let's show the results for first bin
plot_im2D(abs.(x_approx[:,:,1,1,1,1]),title = "Reco from AcquisitionData")


# We can extract the cartesian kspace (works for 2D or 3D) with the function `kDataCart`
kspace = kDataCart(acqData)
size(kspace)

# Dimensions of the kspace are : kx, ky, kz, Channels, Echoes, Repetitions

# Let's see a standard reconstruction with the `ifft` function works :
x_approx2 = ifftshift(ifft(ifftshift(kspace)))
plot_im2D(abs.(x_approx2[:,:,1,1,1,1]),title = "Reco from extracted kspace")

# We can create create and apply a mask on the kspace
mask = ones(eltype(kspace),size(kspace))
mask[1:2:end,:,:,:,:,:] .= 0
acq_u = AcquisitionData(kspace .* mask)
params[:reco] = "direct"
params[:reconSize] = (N,N2)

x_u = reconstruction(acq_u, params)
plot_im2D(abs.(x_u[:,:,1,1,1,1]),title = "Reco from undersampled AcquisitionData")

# ## Reproducibility

# This page was generated with the following version of Julia:

using InteractiveUtils
io = IOBuffer();
versioninfo(io);
split(String(take!(io)), '\n')

# And with the following package versions

import Pkg; Pkg.status()
