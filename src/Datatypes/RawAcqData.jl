export RawAcquisitionData, EncodingCounters, AcquisitionHeader

struct RawAcquisitionData
  params::Dict{String, Any}
  profiles::Vector{Profile}
end

struct Profile
  head::AcquisitionHeader
  traj::Array{Float32,2}
  data::Array{Complex{Float32},2}
end

struct AcquisitionHeader
  version::Int16
  flags::Int64
  measurement_uid::Int32
  scan_counter::Int32
  acquisition_time_stamp::Int32
  physiology_time_stamp::Vector{Int32}
  number_of_samples::Int16
  available_channels::Int16
  active_channels::Int16
  channel_mask::Vector{Int64}
  discard_pre::Int16
  discard_post::Int16
  center_sample::Int16
  encoding_space_ref::Int16
  trajectory_dimensions::Int16
  sample_time_us::Float32
  position::Vector{Float32}
  read_dir::Vector{Float32}
  phase_dir::Vector{Float32}
  slice_dir::Vector{Float32}
  patient_table_position::Vector{Float32}
  idx::EncodingCounters
  user_int::Vector{Int32}
  user_float::Vector{Float32}
end

struct EncodingCounters
  kspace_encode_step_1::Int16
  kspace_encode_step_2::Int16
  average::Int16
  slice::Int16
  contrast::Int16
  phase::Int16
  repetition::Int16
  set::Int16
  segment::Int16
  user::Vector{Int16}
end
