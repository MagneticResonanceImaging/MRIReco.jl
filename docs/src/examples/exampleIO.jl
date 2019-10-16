using HTTP, PyPlot, MRIReco

if !isdir("brukerfileCart")
  HTTP.open("GET", "http://media.tuhh.de/ibi/mrireco/brukerfileCart.zip") do http
    open("brukerfileCart.zip", "w") do file
        write(file, http)
    end
  end
  run(`unzip brukerfileCart.zip`)
end

filename = "simple_spiral.h5"
if !isfile(filename)
  HTTP.open("GET", "https://netix.dl.sourceforge.net/project/ismrmrd/data/simple_spiral.h5") do http
    open(filename, "w") do file
      write(file, http)
    end
  end
end

f = BrukerFile("brukerfileCart")
raw = RawAcquisitionData(f)
acq = AcquisitionData(raw)

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (acq.encodingSize[1],acq.encodingSize[2])

img = reconstruction(acq, params)

filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/bruker.png")
exportImage(filename, abs.(img[:,:,10,1,1]) )


fout = ISMRMRDFile(@__DIR__()*"outputfile.h5")
save(fout, raw)
