#---------------------------------------------------------
# # [Binning](@id 01-binning)
#---------------------------------------------------------


#=
## Description

This example described how to perform data binning with different number of profiles
which is generaly used for self-gating acquisition
in order to reconstruct images along the cardiac/respiratory cycle :
1. Trotier AJ, Castets CR, Lefrançois W, et al. USPIO-enhanced 3D-cine self-gated cardiac MRI based on a stack-of-stars golden angle short echo time sequence: Application on mice with acute myocardial infarction. Journal of Magnetic Resonance Imaging 2016;44:355–365 doi: 10.1002/jmri.25150.
2. Ribot EJ, Duriez TJ, Trotier AJ, Thiaudiere E, Franconi J-M, Miraux S. Self-gated bSSFP sequences to detect iron-labeled cancer cells and/or metastases in vivo in mouse liver at 7 Tesla. Journal of Magnetic Resonance Imaging 2015;41:1413–1421 doi: 10.1002/jmri.24688.


Here, we will create a simulated 2D radial acquisition and split the projections in 2 parts along the contrast dimension.
The number of projection into each bin will be different to show how MRIReco handle that case.
=#

# ## Setup
using CairoMakie
using ImageUtils: shepp_logan
using MRIReco, MRISimulation

function plot_im2D(im2D;title::String="")
    f = Figure()
    ax = Axis(f[1, 1],aspect = DataAspect(), yreversed = true, title = title)
    image!(ax, im2D')
    hidedecorations!(ax, grid = false)
    f
end

# ## Simulate a radial acquisition
N = 256
x = shepp_logan(N)

params = Dict{Symbol, Any}()
params[:simulation] = "fast"
params[:trajName] = "Radial"
params[:numProfiles] = round(Int64,400)
params[:numSamplingPerProfile] = round(Int64,N)

acqRad = simulation(x, params)
rawRad = RawAcquisitionData(acqRad)

# For real acquisition we first create the RawAcquisitionData structure
# and then convert into the AcquisitionData.

# ## Binning Data
rawRad2 = deepcopy(rawRad)
for i in 1:length(rawRad.profiles)
    if  mod(i,4) == 0
        rawRad2.profiles[i].head.idx.contrast = 0
    else
        rawRad2.profiles[i].head.idx.contrast = 1
    end
end;

# We need to tell to MRIReco that our trajectory is a custom one :
rawRad2.params["trajectory"] = "custom";

# ## Reconstruction

# To perform the reconstruction we need to convert the RawAcquisitionData into and AcquisitionData structure.

acqRad2 = AcquisitionData(rawRad2)
# We can also plot the number of projection for both bin :
nPro1 = Int32(length(acqRad2.kdata[1,1,1])/N)
nPro2 = Int32(length(acqRad2.kdata[2,1,1])/N)
println("Number of projection in : \n
- Bin 1 = $nPro1\n
- Bin 2 = $nPro2")

# Now we can perform a standard reconstruction
params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N,N)

Ireco = reconstruction(acqRad2, params)
size(Ireco)

# let's show the results for first bin
plot_im2D(abs.(Ireco[:,:,1,1]),title = "First bin")

# and the second bin
plot_im2D(abs.(Ireco[:,:,1,2]),title = "Second bin")

# As expected we have more streaking artifacts on the first bin because we reconstruct the image with less projections.

# ## Reproducibility

# This page was generated with the following version of Julia:

io = IOBuffer();
versioninfo(io);
split(String(take!(io)), '\n')

# And with the following package versions

import Pkg; Pkg.status()
