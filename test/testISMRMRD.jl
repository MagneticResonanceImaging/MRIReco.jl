using MRIReco
using Test
using HTTP

@testset "ISMRMRD" begin

filename = joinpath(datadir, "simple_gre.h5")

f = ISMRMRDFile(filename)
acq = RawAcquisitionData(f)

@test acq.profiles[1].head.version == Int16(0)
@test acq.profiles[1].head.measurement_uid == Int32(37)
@test acq.profiles[1].head.scan_counter == Int32(1)
@test acq.profiles[1].head.number_of_samples == Int16(256)
@test acq.profiles[1].head.available_channels == Int16(32)
@test acq.profiles[1].head.active_channels == Int16(32)

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (256,128) #this should be clear from context

Ireco = reconstruction(AcquisitionData(acq), params)
exportImage(joinpath(tmpdir,"recogre.png"), abs.(Ireco[:,:,1,1,1]))

filenameCopy = joinpath(tmpdir, "simple_gre_copy.h5")
fCopy = ISMRMRDFile(filenameCopy)
# store data in another ISMRMRD file
save(fCopy, acq)
acqCopy = RawAcquisitionData(fCopy)

@test convert(MRIReco.AcquisitionHeaderImmutable, acq.profiles[1].head) == 
       convert(MRIReco.AcquisitionHeaderImmutable, acqCopy.profiles[1].head) 
@test acqCopy.profiles[1].data == acq.profiles[1].data

IrecoCopy = reconstruction(AcquisitionData(acqCopy), params)

@test IrecoCopy == Ireco

filename = joinpath(datadir, "simple_spiral.h5")

f = ISMRMRDFile(filename)

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (128,128) #this should be clear from context

Ireco = reconstruction(AcquisitionData(f), params)
exportImage(joinpath(tmpdir,"recospiral.png"), abs.(Ireco[:,:,1,1,1]))

filenameCopy = joinpath(tmpdir, "simple_spiral_copy.h5")
fCopy = ISMRMRDFile(filenameCopy)
# store data in another ISMRMRD file
acq = RawAcquisitionData(f)
save(fCopy, acq)
acqCopy = RawAcquisitionData(f)

@test convert(MRIReco.AcquisitionHeaderImmutable, acq.profiles[1].head) == 
       convert(MRIReco.AcquisitionHeaderImmutable, acqCopy.profiles[1].head) 
@test acqCopy.profiles[1].traj == acq.profiles[1].traj
@test acqCopy.profiles[1].data == acq.profiles[1].data

# 3D_partial_fourier.h5

#=filename = "3D_partial_fourier.h5"
if !isfile(filename)
  HTTP.open("GET", "https://netix.dl.sourceforge.net/project/ismrmrd/data/3D_partial_fourier.h5") do http
    open(filename, "w") do file
      write(file, http)
    end
  end
end

f = ISMRMRDFile(filename)
=#

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
