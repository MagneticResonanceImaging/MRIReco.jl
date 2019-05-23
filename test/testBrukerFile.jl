@testset "BrukerFile" begin

if !isdir("brukerfileCart")
  HTTP.open("GET", "http://media.tuhh.de/ibi/mrireco/brukerfileCart.zip") do http
    open("brukerfileCart.zip", "w") do file
        write(file, http)
    end
  end
  run(`unzip brukerfileCart.zip`)
end

b = BrukerFile("brukerfileCart")

acq = RawAcquisitionData(b)
acqData = AcquisitionData(acq)
N = acqData.encodingSize

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N[1],N[2]) #this should be clear from context

Ireco = reconstruction(acqData, params)
exportImage("brukerCart.png", Ireco[:,:,1,1,1])

# Convert to ISMRMRD file
fout = ISMRMRDFile("brukerfileCart.h5")
save(fout, acq)

# Reco data stored in BrukerFile
Iloaded = recoData(b)
@test size(Iloaded) == (128, 128, 15)




end
