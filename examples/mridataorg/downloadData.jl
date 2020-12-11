using HTTP

function download_(filenameServer, filenameLocal)
  if !isfile(filenameLocal)
    @info "download $(filenameLocal)..."
    HTTP.open("GET", filenameServer) do http
      open(filenameLocal, "w") do file
        write(file, http)
      end
    end
  end
end

mkpath("./data/")

# Download data
download_("http://mridata.org/download/fd28ac5d-ed5a-462f-9064-0a5a13d76291", "./data/ksp_knee_3dfse.h5")