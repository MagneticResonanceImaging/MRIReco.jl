export ISMRMRD

# This is one of the most complicated file formats I have seen.

struct ISMRMRD <: MRIFile
  filename::String
  header
end

struct ISMRMRDEncodingCounters
  kspace_encode_step_1::UInt16
  kspace_encode_step_2::UInt16
  average::UInt16
  slice::UInt16
  contrast::UInt16
  phase::UInt16
  repetition::UInt16
  set::UInt16
  segment::UInt16
  user::NTuple{8,UInt16}
end

struct ISMRMRDAcquisitionHeader
  version::UInt16
  flags::UInt64
  measurement_uid::UInt32
  scan_counter::UInt32
  acquisition_time_stamp::UInt32
  physiology_time_stamp::NTuple{3,UInt32}
  number_of_samples::UInt16
  available_channels::UInt16
  active_channels::UInt16
  channel_mask::NTuple{16,UInt64}
  discard_pre::UInt16
  discard_post::UInt16
  center_sample::UInt16
  encoding_space_ref::UInt16
  trajectory_dimensions::UInt16
  sample_time_us::Float32
  position::NTuple{3,Float32}
  read_dir::NTuple{3,Float32}
  phase_dir::NTuple{3,Float32}
  slice_dir::NTuple{3,Float32}
  patient_table_position::NTuple{3,Float32}
  idx::ISMRMRDEncodingCounters
  user_int::NTuple{8,Int32}
  user_float::NTuple{8,Float32}
end

#function Base.write(io::IO, p::ParamsType)
#  write(io, reinterpret(Int8,[p]))
#end

function ISMRMRD(filename::String)
  headerStr = h5read(filename, "/dataset/xml")
  xdoc = parse_string(headerStr[1])

  # This is a compouned type and HD5.jl cannot
  # handle it
  data = h5read(filename, "/dataset/data")
end
