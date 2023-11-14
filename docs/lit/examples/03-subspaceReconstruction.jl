#---------------------------------------------------------
# # [Subspace](@id 03-subspaceReconstruction)
#---------------------------------------------------------

# ## Description
# 
# This example described how to perform a subspace reconstruction for $T_2$ mapping acceleration.

# ## Setup and define global variable
using CairoMakie
using ImageUtils: shepp_logan
using LinearAlgebra
using  Random
using MRIReco, MRISimulation, MRICoilSensitivities, MRISampling,MRIOperators
using RegularizedLeastSquares
color=Makie.wong_colors() # color for plots

N = 128
T = ComplexF32
nCh = 4
nEchos = 10
TE = 7.0



x = T.(shepp_logan(N))

# ## simulate MESE acquisition of a shepp logan phantom
#
# First we need to simulate a multi-echo spin-echo phantom. The R₂ maps is based on the
# value of each voxel the shepp logan phantom. 
# For the T2 we use
rmap = 0.05*abs.(x)
TEnum = Float64.(collect(TE:TE:TE*nEchos))

coilsens = T.(birdcageSensitivity(N, nCh, 4.))
params = Dict{Symbol,Any}()
params[:simulation] = "fast"
params[:trajName] = "Cartesian"
params[:numProfiles] = floor(Int64, N)
params[:numSamplingPerProfile] = N
params[:r2map] = rmap
params[:T_echo] = TEnum
params[:seqName] = "ME"
params[:refocusingAngles] = Float64.(repeat([pi], length(TEnum)))
params[:senseMaps] = coilsens

acqData = simulation(real(x), params)

# ## Subsampling the kspace
# We need to subsample the data. We generate a sampling mask which is different for each TE.

mask = zeros(Int32,(N,N,length(TEnum)))
for echo = 1:length(TEnum)
  mask_tmp = zeros(Int32,(N,N))
  cal_disk = sample_poissondisk((N,N),4.0;calsize = 14,seed=echo)

  mask_tmp[cal_disk] .= 1
  mask[:,:,echo] = mask_tmp
end

f=Figure()
ax = Axis(f[1,1],title = "mask TE n°1")
heatmap!(ax,mask[:,:,1])
ax = Axis(f[1,2],title = "mask TE n°2")
heatmap!(ax,mask[:,:,2])
f

# Let's subsample the k-space in 2 different ways : 
# - With the same mask along the temporal dimensions : `acqData_u`
# - With different masks along the temporal dimensions : `acqData_u2`

mask_multi2 = repeat(mask[:,:,1],1,1,1,nCh,nEchos,1);# same first mask along TE

kspace = kDataCart(acqData);
kspace = kspace .* mask_multi2;
acqData_u = AcquisitionData(kspace);

mask_multi = repeat(mask,1,1,1,nCh,1,1);
mask_multi = permutedims(mask_multi,(1,2,5,4,3,6));

kspace = kDataCart(acqData);
kspace = kspace .* mask_multi;
acqData_u2 = AcquisitionData(kspace);

# ## Reconstruction of undersampled rawdata
params = Dict{Symbol,Any}()
params[:reconSize] = (N, N)
params[:reco] = "multiCoilMultiEcho"
params[:reg] = L2Regularization(1.e-3)
params[:iterations] = 1
params[:solver] = CGNR
params[:senseMaps] = reshape(coilsens, N, N, 1, nCh)

im_x = reconstruction(acqData, params).data # fully reconstruction
im_phant = abs.(im_x)

im_x_u = reconstruction(acqData_u, params).data # undersampling same mask
im_phant_u = abs.(im_x_u)

im_x_u2 = reconstruction(acqData_u2, params).data # undersampling different mask along TE
im_phant_u2 = abs.(im_x_u2)

p1 = (100, 45)
p2 = (40, 60)
echo=1
f = Figure()
ax = Axis(f[1, 1], aspect=1, title="Fully")
hidedecorations!(ax)
heatmap!(ax, im_phant[:, :, 1, echo, 1, 1])
scatter!(ax, p1,color=color[4])
scatter!(ax, p2,color=color[4])

ax = Axis(f[2, 1], aspect=1, title="CS same mask")
hidedecorations!(ax)
heatmap!(ax, im_phant_u[:, :, 1, echo, 1, 1])

ax = Axis(f[3, 1], aspect=1, title="CS different masks")
hidedecorations!(ax)
heatmap!(ax, im_phant_u2[:, :, 1, echo, 1, 1])

ax = Axis(f[:, 2],title = "T2₁=$(round(1/rmap[p1...],digits=0)) | T2₂=$(round(1/rmap[p2...],digits=0)) ms",xlabel = "TE [ms]",ylabel = "MR signal")
lines!(ax, im_phant[p1..., 1, :, 1, 1],color=:black,label = "fully")
lines!(ax, im_phant[p2..., 1, :, 1, 1],color=:black)

lines!(ax, im_phant_u[p1..., 1, :, 1, 1],color=color[2],label = "same mask")
lines!(ax, im_phant_u[p2..., 1, :, 1, 1],color=color[2],)

lines!(ax, im_phant_u2[p1..., 1, :, 1, 1],color=color[3],label = "different mask")
lines!(ax, im_phant_u2[p2..., 1, :, 1, 1],color=color[3])
axislegend(ax)
f



# When the same mask is used along the temporal dimension, we see a standard decaying
# exponential curve. However the rate of decrease is biased by the PSF effect of the
# undersampling mask, corresponding by the same sum of others weigthed pixels. 
#
# When the mask is not the same along the temporal dimension, we observed a noisy curve
# close to the real exponential.



# # Implementation of a subspace reconstruction
# ## Build the dictionary
# We build a signal dictionary using the analytical equation :
#
# $$S(TE) = \frac{-TE}{T2}$$ 
#
# with a range of T2 from 1:1:2000 ms

function createExpBasis(TE_vec::AbstractVector{T}, T2_vec::AbstractVector{T}) where {T<:Real}
  nTE = length(TE_vec)
  nSimu = length(T2_vec)
  expSignal = zeros(T, nSimu, nTE)

  for (i, T2) in enumerate(T2_vec)
    expSignal[i, :] = exp.(-TE_vec / T2_vec[i])
  end

  return expSignal
end

T2_vec = (1:1:2000)
sDict = createExpBasis(Float32.(TEnum), Float32.(T2_vec))

# For real application with stimulated echo we have to remove the 1st echo from the dictionary and the rawdata

# ## Extract the subspace temporal basis
# Now we can perform an svd decomposition of the signal dictionnary in order to extract the temporal basis in order to use them during the reconstruction
svd_obj = svd(sDict)
basis = Complex.(svd_obj.V)[:, 1:5]
# In this example, we will use the 5 first basis. 
f=Figure()
ax = Axis(f[1,1])
for i = 1:size(basis,2)
  lines!(ax,real.(svd_obj.V[:,i]))
end
f
# Our supposition for the subspace reconstruction is that the MESE signal can be
# approximated by a linear combinaison of the 5 basis corresponding to the temporal curve.

# ## Subspace reconstruction
# In order to perform the subspace reconstruction we are using a dedicated pipeline name :
# `params[:reco] = "multiCoilMultiEchoSubspace"` and we also need to pass the basis
# `params[:basis] = basis`
# 
# Here, we also have added a Wavelet spatial regularization that will be applied on the basis coefficient maps.

params = Dict{Symbol,Any}()
params[:reconSize] = (N, N)
params[:reco] = "multiCoilMultiEchoSubspace"

params[:reg] = L1Regularization(0.001)
params[:sparseTrafo] = "Wavelet" #sparse trafo
params[:solver] = ADMM
params[:senseMaps] = reshape(coilsens, N, N, 1, nCh)
params[:basis] = basis

params[:iterations] = 1
α = reconstruction(acqData, params)
im_sub = abs.(applySubspace(α, basis))

params[:iterations] = 100
α = reconstruction(acqData_u, params)
im_sub_2 = abs.(applySubspace(α, basis))

α = reconstruction(acqData_u2, params)
im_sub_3 = abs.(applySubspace(α, basis));

# ## Coefficient maps
# The reconstruction returns the coefficient maps :

f = Figure()
for i = 1:size(basis,2)
  ax = Axis(f[1,i],aspect=1,title = "Basis n°$i")
  heatmap!(abs.(α[:,:,1,i,1,1]),colormap=:plasma)
  hidedecorations!(ax)
end
f

# ## Virtual echo images
# We need to multiply the subspace basis to the coefficient maps in order to get the virtual TE images
#
# $$im_{TE}(i,j) = \sum_{basis=1}^{5} \alpha(i,j) \times basis'$$
#
# which gives the following results :
begin
p1 = (100, 45)
p2 = (40, 60)
echo = 1
f = Figure()
ax = Axis(f[1, 1], aspect=1, title="Fully standard")
hidedecorations!(ax)
heatmap!(ax, im_phant[:, :, 1, echo, 1, 1])
scatter!(ax, p1,color=:red)
scatter!(ax, p2,color=:red)

ax = Axis(f[1, 2], aspect=1, title="Fully with subspace reco")
hidedecorations!(ax)
heatmap!(ax, im_sub[:, :, 1, echo, 1, 1])

ax = Axis(f[1, 3], aspect=1, title="CS same mask")
hidedecorations!(ax)
heatmap!(ax, im_sub_2[:, :, 1, echo, 1, 1])

ax = Axis(f[1, 4], aspect=1, title="CS different masks")
hidedecorations!(ax)
heatmap!(ax, im_sub_3[:, :, 1, 1, 1, 1])

ax = Axis(f[2, :],title = "T2₁=$(round(1/rmap[p1...],digits=0)) | T2₂=$(round(1/rmap[p2...],digits=0)) ms",xlabel = "TE [ms]",ylabel = "MR signal")
function scale_T2(T2vect)
  return T2vect/T2vect[1]
end
lines!(ax, scale_T2(im_phant[p1..., 1, :, 1, 1]),color=:black,label = "fully")
lines!(ax, scale_T2(im_phant[p2..., 1, :, 1, 1]),color=:black)

lines!(ax, scale_T2(im_sub[p1..., 1, :, 1, 1]),color=color[1],label = "fully subspace")
lines!(ax, scale_T2(im_sub[p2..., 1, :, 1, 1]),color=color[1])

lines!(ax, scale_T2(im_sub_2[p1..., 1, :, 1, 1]),color=color[2],label = "same mask")
lines!(ax, scale_T2(im_sub_2[p2..., 1, :, 1, 1]),color=color[2])

lines!(ax, scale_T2(im_sub_3[p1..., 1, :, 1, 1]),color=color[3],label = "different mask")
lines!(ax, scale_T2(im_sub_3[p2..., 1, :, 1, 1]),color=color[3])
axislegend(ax)
f
end
# # Reproducibility

# This page was generated with the following version of Julia:

using InteractiveUtils
io = IOBuffer();
versioninfo(io);
split(String(take!(io)), '\n')

# And with the following package versions

import Pkg; Pkg.status()