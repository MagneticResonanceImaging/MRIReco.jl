using CairoMakie, MRIReco, MRIReco.RegularizedLeastSquares
using MRIFiles, MRISampling, MRICoilSensitivities
include(joinpath(@__DIR__,"exampleUtils.jl"))
# load fully sampled data
f = ISMRMRDFile(@__DIR__()*"/data/knee_3dFSE_slice170.h5")
acqData = AcquisitionData(f);

# reconstruct
params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (320,320) # this size is also contained in acqData.encodingSize

img = reconstruction(acqData, params)
img = sqrt.(sum(img.^2,dims=5))

# To simulate an undersampled reconstruction, we retrospectively undersample the data using a Poisson disk pattern.

redFac = 4.0
acqDataSub = sample_kspace(acqData,redFac,"poisson",calsize=30,profiles=false);

# show sampling pattern
msk = zeros(acqDataSub.encodingSize[1],acqDataSub.encodingSize[2])
msk[acqDataSub.subsampleIndices[1]] .= 1

# Estimate the coil sensitivities using ESPIRiT and reconstruct using SENSE

# coil sensitivities
smaps = espirit(acqData,(6,6),30,eigThresh_2=0)

# SENSE reconstruction
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (320,320)
params[:senseMaps] = smaps

params[:solver] = CGNR
params[:reg] = L2Regularization(1.e-4)
params[:iterations] = 5
params[:normalizeReg] = MeasurementBasedNormalization()

img_cg = reconstruction(acqDataSub, params)

# Using TV regularization recquires us to change some parameters.

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (320,320)
params[:senseMaps] = smaps

params[:solver] = ADMM
params[:reg] = TVRegularization(2e-1, shape = (320, 320))
params[:iterations] = 50
params[:rho] = 0.1
params[:absTol] = 1.e-4
params[:relTol] = 1.e-2
params[:tolInner] = 1.e-2
params[:normalizeReg] = MeasurementBasedNormalization()


img_tv = reconstruction(acqDataSub, params)


begin
  colormap=:grays
  f = Figure(size=(800,800))
  ax = Axis(f[1,1],title="Phantom")
  heatmap!(ax,rotr90(abs.(img[:,:,1,1,1]));colormap)
  ax = Axis(f[1,2],title="Mask")
  heatmap!(ax,rotr90(abs.(msk));colormap)
  ax = Axis(f[2,1],title="CG Reconstruction")
  heatmap!(ax,rotr90(abs.(img_cg[:,:,1,1,1]));colormap)
  ax = Axis(f[2,2],title="TV Reconstruction")
  heatmap!(ax,rotr90(abs.(img_tv[:,:,1,1,1]));colormap)
  [hidedecorations!(f.content[ax]) for ax in eachindex(f.content)]
  f
end


# export images
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/kneeOrig.png")
exportImage(filename, abs.(img[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/mask.png")
exportImage(filename, abs.(msk[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/kneeCG.png")
exportImage(filename, abs.(img_cg[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/kneeTV.png")
exportImage(filename, abs.(img_tv[:,:,1,1,1]) )