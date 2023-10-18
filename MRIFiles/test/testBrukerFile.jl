@testset "BrukerFile" begin

@testset "BrukerFile read parameters" begin
    b = BrukerFile( joinpath(datadir, "BrukerFile", "2D_RARE") )
    t1 = @elapsed b["ExcPulse1"]

    @test b["ExcPulse1"] == "(1.3125, 3200, 90, Yes, 3, 4200, 0.236151639875348, 0.200434548747244, 0, 50, 0.317196887605166, <\$ExcPulse1Shape>)"
    @test b["PVM_Fov"] == ["27","27"]
    @test b["PVM_FreqDriftYN"] == "Yes"

    @test b["VisuCoreWordType"] == "_16BIT_SGN_INT"
    @test b["VisuFGOrderDesc"][1] == Any[15.0, "<FG_SLICE>", "<>", 0.0, 2.0]
    @test b["VisuCoreDataOffs"] == repeat(["0","0","0"],5)

    t2 = @elapsed b["ExcPulse1"]
    @info "Read parameters timing" t1 t2
    @test t2 < 0.0001
end

@testset "BrukerFile read 2dseq" begin
    b = BrukerFile( joinpath(datadir, "BrukerFile", "2D_FLASH") )

    ## standard reco : Int16 + Magnitude
    ima1 = recoData(b,1)
    ima1 = ima1 ./ maximum(ima1);

    ## Float32 + Magnitude
    ima4 = recoData(b,4)
    ima4 = ima4 ./ maximum(ima4);
    @test MRIReco.norm(vec(ima1)-vec(ima4))/MRIReco.norm(vec(ima1)) < 0.02

    ## complex + Int16 + shuffle Reco
    ima6 = recoData(b,6)
    ima6 = complex.(ima6[:,:,1:4],ima6[:,:,5:end])
    ima6 = sqrt.(sum(abs.(ima6).^2,dims=3));
    ima6 = ima6 ./ maximum(ima6);
    @test MRIReco.norm(vec(ima1)-vec(ima6))/MRIReco.norm(vec(ima1)) < 0.02

    ## UInt8 Magnitude
    ima7 = recoData(b,7)
    ima7 = ima7 ./ maximum(ima7);
    @test MRIReco.norm(vec(ima1)-vec(ima7))/MRIReco.norm(vec(ima1)) < 0.12

end

@testset "BrukerFile Reco" begin
b = BrukerFile( joinpath(datadir, "BrukerFile", "2D_RARE") )

acq = RawAcquisitionData(b)
acqData = AcquisitionData(acq)

params = Dict{Symbol, Any}()
params[:reco] = "direct"

Ireco = reconstruction(acqData, params)
exportImage( joinpath(tmpdir, "brukerCart.png"), abs.(Ireco[:,:,1,1,1]))

# Convert to ISMRMRD file
fout = ISMRMRDFile(joinpath(tmpdir, "brukerfileCart.h5"))
save(fout, acq)

# Reco data stored in BrukerFile
Iloaded = recoData(b)
@test size(Iloaded) == (128, 128, 15)

## Test reconstruction for multi-coil datasets (2D and 3D FLASH)
listBrukFiles = ["2D_FLASH", "3D_FLASH"]
listNormValues = [0.02, 0.15]

for i = 1:length(listBrukFiles)
    b = BrukerFile( joinpath(datadir, "BrukerFile", listBrukFiles[i]) )
    raw = RawAcquisitionData(b)
    acq = AcquisitionData(raw)
    params = Dict{Symbol, Any}()
    params[:reco] = "direct"
    
    Ireco = reconstruction(acq, params);
    @test size(Ireco) == (raw.params["encodedSize"][1], raw.params["encodedSize"][2], raw.params["encodedSize"][3], numContrasts(acq), raw.params["receiverChannels"], numRepetitions(acq))

    Isos = sqrt.(sum(abs.(Ireco).^2,dims=5));
    Isos = Isos ./ maximum(Isos);

    I2dseq = recoData(b)
    I2dseq = I2dseq ./ maximum(I2dseq);

    @test norm(vec(I2dseq)-vec(Isos))/norm(vec(I2dseq)) < listNormValues[i]
end


## Test reconstruction for multi-coil datasets with CS acceleration
b = BrukerFile( joinpath(datadir, "BrukerFile", "PV360_3D_FLASH_CS" ))
raw = RawAcquisitionData(b)
acq = AcquisitionData(raw,OffsetBruker=true)

sens = espirit(acq,eigThresh_2 = 0);

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"    # encoding model
params[:reconSize] = acq.encodingSize
params[:sparseTrafo] = "Wavelet" #sparse trafo
params[:regularization] = "L1"       # regularization
params[:λ] = 30.0
params[:solver] = "admm"    # solver
params[:iterations] = 3
params[:iterationsInner] = 2
params[:ρ] = 1.0e-1
params[:senseMaps] = sens

img_CS = abs.(reconstruction(acq, params))
img_CS = permutedims(img_CS[:,:,:,1,1,1],[2,1,3])
img_CS = img_CS ./ maximum(img_CS);

I2dseq = recoData(b)
I2dseq = I2dseq ./ maximum(I2dseq);

@test norm(vec(I2dseq)-vec(img_CS))/norm(vec(I2dseq)) < 0.2


if !Sys.iswindows()
# Reconstruction of 3DUTE
@info "Reconstruction of 3DUTE"
b = BrukerFile( joinpath(datadir, "BrukerFile", "3D_UTE_NR2") )

raw = RawAcquisitionData(b);
acq = AcquisitionData(raw);  # TODO vérification des modifications qui ont étaient effectuée

params = Dict{Symbol, Any}()
params[:reco] = "direct"

Ireco = reconstruction(acq, params);

Isos = sqrt.(sum(abs.(Ireco).^2,dims=5));
Isos = Isos ./ maximum(Isos);
I2dseq = recoData(b)
I2dseq = I2dseq ./ maximum(I2dseq);
# reorient
I2dseq = permutedims(I2dseq,(2,1,3,4))
I2dseq = circshift(I2dseq,(0,0,1,0))

#test reconstruction of repetitions
@test size(Isos,6)==2
@test MRIReco.norm(vec(I2dseq)-vec(Isos[:,:,:,1,1,:]))/MRIReco.norm(vec(I2dseq)) < 0.1

end

end

end
