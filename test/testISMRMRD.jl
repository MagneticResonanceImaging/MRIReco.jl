using HTTP

@testset "ISMRMRD" begin

filename = "simple_gre.h5"

if !isfile(filename)
  HTTP.open("GET", "https://netix.dl.sourceforge.net/project/ismrmrd/data/simple_gre.h5") do http
    open(filename, "w") do file
      write(file, http)
    end
  end
end

f = ISMRMRD(filename)

@test f.head[1].version == Int16(0)
@test f.head[1].measurement_uid == Int32(37)
@test f.head[1].scan_counter == Int32(1)
@test f.head[1].number_of_samples == Int16(256)
@test f.head[1].available_channels == Int16(32)
@test f.head[1].active_channels == Int16(32)

@test size(f.data) == (256, 32, 1281)

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:shape] = (256,128) #this should be clear from context

Ireco = abs.(vec(reconstruction(acquisitionData(f), params)))
Icolored = colorview(Gray, Ireco./maximum(Ireco))
save("recogre.png", reshape(Icolored,256,128,32)[:,:,1] )


filename = "simple_spiral.h5"
if !isfile(filename)
  HTTP.open("GET", "https://netix.dl.sourceforge.net/project/ismrmrd/data/simple_spiral.h5") do http
    open(filename, "w") do file
      write(file, http)
    end
  end
end

f = ISMRMRD(filename)
@test size(f.data) == (1000, 32, 160)



# 3D_partial_fourier.h5

end

#=

For our own reference: Here is code that can be used with Python to load an
ISMRMRD file

import h5py
f = h5py.File('simple_gre.h5', 'r')
dset = f['/dataset/data']
dset[0]['head']

=#
