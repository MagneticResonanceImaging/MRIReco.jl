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

## Test reconstruction for multi-coil datasets (2D and 3D FLASH)
listBrukFiles = ["2D_FLASH","3D_FLASH"]
listNormValues = [0.02, 0.15]

for i = 1:length(listBrukFiles)
    b = BrukerFile( joinpath(datadir, "BrukerFile", listBrukFiles[i]) )
    raw = RawAcquisitionData(b)
    acq = AcquisitionData(raw)
    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    if (acq.encodingSize[3]>1)
        params[:reconSize] = (acq.encodingSize[1],acq.encodingSize[2],acq.encodingSize[3]);
    else
        params[:reconSize] = (acq.encodingSize[1],acq.encodingSize[2]);
    end
    Ireco = reconstruction(acq, params);
    @test size(Ireco) == (raw.params["encodedSize"][1], raw.params["encodedSize"][2], raw.params["encodedSize"][3], 1, raw.params["receiverChannels"])

    Isos = sqrt.(sum(abs.(Ireco).^2,dims=5));
    Isos = Isos ./ maximum(Isos);

    I2dseq = recoData(b)
    I2dseq = I2dseq ./ maximum(I2dseq);

    @test norm(vec(I2dseq)-vec(Isos))/norm(vec(I2dseq)) < listNormValues[i]
end
end
