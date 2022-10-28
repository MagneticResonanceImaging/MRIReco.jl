using PyPlot, MRIReco

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
smaps = espirit(acqData,(6,6),30,eigThresh_1=0.035,eigThresh_2=0.98)

# SENSE reconstruction
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (320,320)
params[:senseMaps] = smaps

params[:solver] = "cgnr"
params[:regularization] = "L2"
params[:λ] = 1.e-4
params[:iterations] = 5
params[:normalizeReg] = true

img_cg = reconstruction(acqDataSub, params)

# Using TV regularization recquires us to change some parameters.

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (320,320)
params[:senseMaps] = smaps

params[:solver] = "admm"
params[:regularization] = "TV"
params[:λ] = 1.e-1 # 5.e-2
params[:iterations] = 50
params[:ρ] = 0.1
params[:absTol] = 1.e-4
params[:relTol] = 1.e-2
params[:tolInner] = 1.e-2
params[:normalizeReg] = true

img_tv = reconstruction(acqDataSub, params)


# use PyPlot for interactive display
figure(1)
clf()
subplot(2,2,1)
title("Phantom")
imshow(abs.(img[:,:,1,1,1]))
subplot(2,2,2)
title("Mask")
imshow(abs.(msk))
subplot(2,2,3)
title("CG Reconstruction")
imshow(abs.(img_cg[:,:,1,1,1]))
subplot(2,2,4)
title("TV Reconstruction")
imshow(abs.(img_tv[:,:,1,1,1]))

# export images
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/kneeOrig.png")
exportImage(filename, abs.(img[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/mask.png")
exportImage(filename, abs.(msk[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/kneeCG.png")
exportImage(filename, abs.(img_cg[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/kneeTV.png")
exportImage(filename, abs.(img_tv[:,:,1,1,1]) )
