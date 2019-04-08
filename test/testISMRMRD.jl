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
@test size(f.data) == (909, 32, 160)

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


end

#=

For our own reference: Here is code that can be used with Python to load an
ISMRMRD file

import h5py
f = h5py.File('simple_gre.h5', 'r')
dset = f['/dataset/data']
dset[0]['head']

=#
