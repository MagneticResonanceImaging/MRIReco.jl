@testset "BrukerFile" begin


b = BrukerFile( joinpath(datadir, "BrukerFile", "2D_RARE") )

acq = RawAcquisitionData(b)
acqData = AcquisitionData(acq)
N = acqData.encodingSize

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N[1],N[2]) #this should be clear from context

Ireco = reconstruction(acqData, params)
exportImage( joinpath(tmpdir, "brukerCart.png"), abs.(Ireco[:,:,1,1,1]))

# Convert to ISMRMRD file
fout = ISMRMRDFile(joinpath(tmpdir, "brukerfileCart.h5"))
save(fout, acq)

# Reco data stored in BrukerFile
Iloaded = recoData(b)
@test size(Iloaded) == (128, 128, 15)




end
