export ISMRMRDFile, saveasISMRMRDFile

include("HDF5Types.jl")

struct ISMRMRDFile
  filename::String
end

function RawAcquisitionData(f::ISMRMRDFile, dataset="dataset")

  h5open(f.filename) do h
    headerStr = read(h["/$(dataset)/xml"])
    params = GeneralParameters(headerStr[1])

    d = read(h["/$(dataset)/data"])
    dtypehead = h5t_get_member_type(datatype(h["/$(dataset)/data"]),0)

    M = length(d)

    profiles = Profile[]

    for m=1:M
      head = read_header(d[m].data[1], dtypehead)
      D = Int(head.trajectory_dimensions)
      chan = Int(head.active_channels)

      N = reinterpret(Int64, d[m].data[2][1:8])[1]
      if N > 0
        ptr = reinterpret(Int64, d[m].data[2][9:end])[1]
        traj = unsafe_wrap(Array,Ptr{Float32}(ptr),(D,div(N,D)))
      else
        traj = Matrix{Float32}(undef,0,0)
      end

      N = reinterpret(Int64, d[m].data[3][1:8])[1]
      if N > 0
        ptr = reinterpret(Int64, d[m].data[3][9:end])[1]
        U = unsafe_wrap(Array,Ptr{ComplexF32}(ptr),(div(N,2),))
        dat = reshape(U,div(N,2*chan),chan)
      else
        dat = Matrix{ComplexF32}(undef,0,0)
      end

      push!(profiles, Profile(head,traj,dat) )
    end
    return RawAcquisitionData(params, profiles)
  end
end

function AcquisitionData(f::ISMRMRDFile, dataset="dataset")
  return AcquisitionData(RawAcquisitionData(f,dataset))
end

function read_header(header, dtypehead)
  # very important to use the offset here, since the offset varies from
  # file to file

  buf = IOBuffer(header)

  seek(buf, h5t_get_member_offset(dtypehead, 0))
  version = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 1))
  flags = read(buf,Int64)
  seek(buf, h5t_get_member_offset(dtypehead, 2))
  measurement_uid = read(buf,Int32)
  seek(buf, h5t_get_member_offset(dtypehead, 3))
  scan_counter = read(buf,Int32)
  seek(buf, h5t_get_member_offset(dtypehead, 4))
  acquisition_time_stamp = read(buf,Int32)
  seek(buf, h5t_get_member_offset(dtypehead, 5))
  physiology_time_stamp = ntuple(i->read(buf,Int32),3)
  seek(buf, h5t_get_member_offset(dtypehead, 6))
  number_of_samples = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 7))
  available_channels = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 8))
  active_channels = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 9))
  channel_mask = ntuple(i->read(buf,Int64),16)
  seek(buf, h5t_get_member_offset(dtypehead, 10))
  discard_pre = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 11))
  discard_post = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 12))
  center_sample = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 13))
  encoding_space_ref = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 14))
  trajectory_dimensions = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtypehead, 15))
  sample_time_us = read(buf,Float32)
  seek(buf, h5t_get_member_offset(dtypehead, 16))
  position = ntuple(i->read(buf,Float32),3)
  seek(buf, h5t_get_member_offset(dtypehead, 17))
  read_dir = ntuple(i->read(buf,Float32),3)
  seek(buf, h5t_get_member_offset(dtypehead, 18))
  phase_dir = ntuple(i->read(buf,Float32),3)
  seek(buf, h5t_get_member_offset(dtypehead, 19))
  slice_dir = ntuple(i->read(buf,Float32),3)
  seek(buf, h5t_get_member_offset(dtypehead, 20))
  patient_table_position = ntuple(i->read(buf,Float32),3)
  ###
  dtype_enc = h5t_get_member_type(dtypehead,21)
  off = h5t_get_member_offset(dtypehead, 21)
  seek(buf, h5t_get_member_offset(dtype_enc, 0) + off)
  kspace_encode_step_1 = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 1) + off)
  kspace_encode_step_2 = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 2) + off)
  average = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 3) + off)
  slice = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 4) + off)
  contrast = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 5) + off)
  phase = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 6) + off)
  repetition = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 7) + off)
  set = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 8) + off)
  segment = read(buf,Int16)
  seek(buf, h5t_get_member_offset(dtype_enc, 9) + off)
  user = ntuple(i->read(buf,Int16),8)
  ###
  idx = EncodingCounters(kspace_encode_step_1,kspace_encode_step_2,
                      average,slice,contrast,phase,repetition,set,segment,user)

  seek(buf, h5t_get_member_offset(dtypehead, 22))
  user_int = ntuple(i->read(buf,Int32),8)
  seek(buf, h5t_get_member_offset(dtypehead, 23))
  user_float = ntuple(i->read(buf,Float32),8)

  head = AcquisitionHeader( version, flags, measurement_uid,
    scan_counter, acquisition_time_stamp, physiology_time_stamp, number_of_samples,
    available_channels, active_channels, channel_mask, discard_pre, discard_post,
    center_sample, encoding_space_ref, trajectory_dimensions,sample_time_us,
    position, read_dir, phase_dir, slice_dir, patient_table_position, idx, user_int,
    user_float
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
