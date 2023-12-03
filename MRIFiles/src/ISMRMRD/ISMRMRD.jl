export ISMRMRDFile, AcquisitionHeaderImmutable, EncodingCountersImmutable

include("HDF5Types.jl")
include("Parameters.jl")

struct ISMRMRDFile <: MRIFile
  filename::String
end

"""
    RawAcquisitionData(f::ISMRMRDFile, dataset="dataset")

reads the `ISMRMRDFile` f and stores the result in a `RawAcquisitionDataObject`
"""
function MRIBase.RawAcquisitionData(f::ISMRMRDFile, dataset="dataset"; 
                                    slice=nothing, repetition=nothing, contrast=nothing)

  h5open(f.filename) do h
    headerStr = read(h["/$(dataset)/xml"])
    params = GeneralParameters(headerStr[1])

    M = size(h["/$(dataset)/data"]) # for some reason this is a matrix 

    # In the next lines we only read the headers of the profiles such that
    # we can later filter on them. We first read the data into a binary blob (NTuple)
    # and the exploit the read implementation that reconstructs the struct field by field
    dTypeHeader = get_hdf5type_acquisition_header_only()
    headerData_ = Array{NTuple{ISMRMRD_HEADER_SIZE,UInt8}}(undef, M)
    HDF5.API.h5d_read(h["/$(dataset)/data"], dTypeHeader, H5S_ALL, H5S_ALL, H5P_DEFAULT, headerData_)
    headerData = [ read(IOBuffer(collect(headerData_[m])),AcquisitionHeader) for m in CartesianIndices(M)]

    profiles = Profile[]

    for m in CartesianIndices(M)
      if (repetition == nothing || headerData[m].idx.repetition in repetition) &&
          (slice == nothing || headerData[m].idx.slice in slice) &&
          (contrast == nothing || headerData[m].idx.contrast in contrast)
         d = h["/$(dataset)/data"][m] # Here we read the header again. This could/should be avoided
         head = read_header(d.head)
         D = Int(head.trajectory_dimensions)
         chan = Int(head.active_channels)
         traj = isempty(d.traj) ? Matrix{Float32}(undef, 0, 0) : reshape(d.traj, D, :)

         if !isempty(d.data)
           dat = reshape(reinterpret(ComplexF32, d.data), :, chan)
         else
           dat = Matrix{ComplexF32}(undef,0,0)
         end

         push!(profiles, Profile(head, traj, dat) )
       end
    end
    return RawAcquisitionData(params, profiles)
  end
end

"""
    AcquisitionData(f::ISMRMRDFile, dataset="dataset")

reads the `ISMRMRDFile` f and stores the result in an `AcquisitionDataObject`
"""
function MRIBase.AcquisitionData(f::ISMRMRDFile, dataset="dataset")
  return AcquisitionData(RawAcquisitionData(f,dataset))
end

function read_header(h)

  idx = EncodingCounters(h.idx.kspace_encode_step_1, h.idx.kspace_encode_step_2,
                      h.idx.average, h.idx.slice, h.idx.contrast, h.idx.phase, h.idx.repetition,
                      h.idx.set, h.idx.segment, ntuple(i->h.idx.user[i],8))

  head = AcquisitionHeader( h.version, h.flags, h.measurement_uid,
    h.scan_counter, h.acquisition_time_stamp,
    ntuple(i->h.physiology_time_stamp[i],3),
    h.number_of_samples, h.available_channels, h.active_channels,
    ntuple(i->h.channel_mask[i],16),
    h.discard_pre, h.discard_post,
    h.center_sample, h.encoding_space_ref, h.trajectory_dimensions, h.sample_time_us,
    ntuple(i->h.position[i],3),
    ntuple(i->h.read_dir[i],3),
    ntuple(i->h.phase_dir[i],3),
    ntuple(i->h.slice_dir[i],3),
    ntuple(i->h.patient_table_position[i],3),
    idx,
    ntuple(i->h.user_int[i],8),
    ntuple(i->h.user_float[i],8)
    )

  return head
end


function FileIO.save(f::ISMRMRDFile, acq::RawAcquisitionData, dataset="dataset")
  h5open(f.filename, "w") do file
    headerStr = GeneralParametersToXML(acq.params)
    write(file, "/$(dataset)/xml", [headerStr]) # [] ensures that we store as var string
    writeProfiles(file, "/$(dataset)/data", acq.profiles)
  end
end 
