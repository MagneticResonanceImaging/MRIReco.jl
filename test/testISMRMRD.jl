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

@test f.profiles[1].head.version == Int16(0)
@test f.profiles[1].head.measurement_uid == Int32(37)
@test f.profiles[1].head.scan_counter == Int32(1)
@test f.profiles[1].head.number_of_samples == Int16(256)
@test f.profiles[1].head.available_channels == Int16(32)
@test f.profiles[1].head.active_channels == Int16(32)

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:shape] = (256,128) #this should be clear from context

Ireco = abs.(reconstruction(acquisitionData(f), params))
Icolored = colorview(Gray, Ireco[:,:,1,1,1]./maximum(Ireco[:,:,1,1,1]))
save("recogre.png", Icolored )


filename = "simple_spiral.h5"
if !isfile(filename)
  HTTP.open("GET", "https://netix.dl.sourceforge.net/project/ismrmrd/data/simple_spiral.h5") do http
    open(filename, "w") do file
      write(file, http)
    end
  end
end

f = ISMRMRD(filename)

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:shape] = (128,128) #this should be clear from context

Ireco = abs.(reconstruction(acquisitionData(f), params))
Icolored = colorview(Gray, Ireco[:,:,1,1,1]./maximum(Ireco[:,:,1,1,1]))
save("recospiral.png", Icolored )


# 3D_partial_fourier.h5

filename = "3D_partial_fourier.h5"
if !isfile(filename)
  HTTP.open("GET", "https://netix.dl.sourceforge.net/project/ismrmrd/data/3D_partial_fourier.h5") do http
    open(filename, "w") do file
      write(file, http)
    end
  end
end

f = ISMRMRD(filename)
#@test size(f.data) == (909, 32, 160)

end

#=

For our own reference: Here is code that can be used with Python to load an
ISMRMRD file

import h5py
f = h5py.File('simple_gre.h5', 'r')
dset = f['/dataset/data']
dset[0]['head']

=#
