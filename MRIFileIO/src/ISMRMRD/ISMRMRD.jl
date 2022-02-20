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
function MRIBase.RawAcquisitionData(f::ISMRMRDFile, dataset="dataset")

  h5open(f.filename) do h
    headerStr = read(h["/$(dataset)/xml"])
    params = GeneralParameters(headerStr[1])

    d = read(h["/$(dataset)/data"])
    M = length(d)
    profiles = Profile[]

    for m=1:M
      head = read_header(d[m].head)
      D = Int(head.trajectory_dimensions)
      chan = Int(head.active_channels)

      traj = isempty(d[m].traj) ? Matrix{Float32}(undef,0,0) : reshape(d[m].traj, D, :)

      if !isempty(d[m].data)
        dat = reshape(reinterpret(ComplexF32,d[m].data), :, chan)
      else
        dat = Matrix{ComplexF32}(undef,0,0)
      end

      push!(profiles, Profile(head,traj,dat) )
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
