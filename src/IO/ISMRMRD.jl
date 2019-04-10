export ISMRMRD

function ISMRMRD(filename::String)
  headerStr = h5read(filename, "/dataset/xml")

  params = GeneralParameters(headerStr[1])

  d = h5read(filename, "/dataset/data")

  M = length(d)

  profiles = Profile[]

  chan = params["receiverChannels"]

  for m=1:M

    head = read_header(d[m].data[1])
    D = Int(head.trajectory_dimensions)

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


function read_header(header)
  buf = IOBuffer(header)

  version = read(buf,Int16)
  flags = read(buf,Int64)
  measurement_uid = read(buf,Int32)
  scan_counter = read(buf,Int32)
  acquisition_time_stamp = read(buf,Int32)
  physiology_time_stamp = read!(buf,zeros(Int32,3))

  # The following line does not really make sense. But there seems to be
  # exactly 20 bytes that need to be read away.
  read!(buf,zeros(Int32,5))

  number_of_samples = read(buf,Int16)
  available_channels = read(buf,Int16)
  active_channels = read(buf,Int16)
  channel_mask = read!(buf,zeros(Int64,16))
  discard_pre = read(buf,Int16)
  discard_post = read(buf,Int16)
  center_sample = read(buf,Int16)
  encoding_space_ref = read(buf,Int16)
  trajectory_dimensions = read(buf,Int16)
  sample_time_us = read(buf,Float32)
  position = read!(buf,zeros(Float32,3))
  read_dir = read!(buf,zeros(Float32,3))
  phase_dir = read!(buf,zeros(Float32,3))
  slice_dir = read!(buf,zeros(Float32,3))
  patient_table_position = read!(buf,zeros(Float32,3))
  ###
  kspace_encode_step_1 = read(buf,Int16)
  kspace_encode_step_2 = read(buf,Int16)
  average = read(buf,Int16)
  slice = read(buf,Int16)
  contrast = read(buf,Int16)
  phase = read(buf,Int16)
  repetition = read(buf,Int16)
  set = read(buf,Int16)
  segment = read(buf,Int16)
  user = read!(buf,zeros(Int16,8))
  ###
  idx = EncodingCounters(kspace_encode_step_1,kspace_encode_step_2,
                      average,slice,contrast,phase,repetition,set,segment,user)

  user_int = read!(buf,zeros(Int32,8))
  user_float = read!(buf,zeros(Float32,8))

  head = AcquisitionHeader( version, flags, measurement_uid,
    scan_counter, acquisition_time_stamp, physiology_time_stamp, number_of_samples,
    available_channels, active_channels, channel_mask, discard_pre, discard_post,
    center_sample, encoding_space_ref, trajectory_dimensions,sample_time_us,
    position, read_dir, phase_dir, slice_dir, patient_table_position, idx, user_int,
    user_float
    )

  return head
end
