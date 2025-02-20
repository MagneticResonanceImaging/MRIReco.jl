using Pkg

isinstalled(pkg::String) = any(x -> x.name == pkg && x.is_direct_dep, values(Pkg.dependencies()))

# Install required packages
for P in ["HTTP", "CairoMakie", "MRIReco", "MRIFiles", "MRICoilSensitivities"]
  !isinstalled(P) && Pkg.add(P)
end

# Download data
include("downloadData.jl")

using CairoMakie, MRIReco, MRIFiles, MRICoilSensitivities, MRIReco.RegularizedLeastSquares, MRIReco.MRIOperators.FFTW

##################
## Data Loading ##
##################
f = ISMRMRDFile("data/ksp_knee_3dfse.h5")
acqRaw = RawAcquisitionData(f)

# kspace-center is not specified in the raw data file.
# So we let MRIReco estimate it
acqData = AcquisitionData(acqRaw; estimateProfileCenter=true)

# convert to 2d
acqData2d = convert3dTo2d(acqData)

# extract slices
sl = [50,100,150,200]
acqData2d.kdata = acqData2d.kdata[:,sl,:]
acqData2d.encodingSize = (274,208)

####################################
## generate coil sensitivity maps ##
####################################
@info "Espirit"
smaps = espirit(acqData2d, (6,6), 30, eigThresh_1=0.04, eigThresh_2=0.98)

#figure(2, figsize=(4,4))
#imshow(abs.(smaps[:,:,4,1]))

###############################
## Reconstruction Parameters ##
###############################
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = Tuple(acqData2d.encodingSize[1:2])
params[:solver] = ADMM
params[:reg] = L1Regularization(2.e-1)
params[:sparseTrafo] = "Wavelet"
params[:iterations] = 30
params[:ρ] = 0.1
params[:absTol] = 1.e-4
params[:relTol] = 1.e-3
params[:tolInner] = 1.e-2
params[:senseMaps] = smaps
params[:normalizeReg] = MeasurementBasedNormalization()

############################
## Perform reconstruction ##
############################

@info "Perform Reconstruction "
@time img = fftshift( reconstruction(acqData2d, params).data[:,:,:,1,1], 2)

###########################
## Visualize the results ##
###########################

f = Figure(size = (190, 600))
for (i, s) in enumerate([1, 2, 3, 4])
  ax = Axis(f[i, 1:4], title = "slice $(sl[s])")
  heatmap!(ax, abs.(img[:, :, s, 1, 1]))
  hidedecorations!(ax)
  hidespines!(ax)
end
f
