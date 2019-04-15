export RawAcquisitionData, EncodingCounters, AcquisitionHeader


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
  user::NTuple{8,Int16}
end

function EncodingCounters(;kspace_encode_step_1=0, kspace_encode_step_2=0,
            average=0, slice=0, contrast=0, phase=0, repetition=0, set=0,
            segment=0, user=ntuple(i->Int16(0),8) )
  return EncodingCounters(kspace_encode_step_1, kspace_encode_step_2,
              average, slice, contrast, phase, repetition, set,
              segment, user)
end

struct AcquisitionHeader
  version::Int16
  flags::Int64
  measurement_uid::Int32
  scan_counter::Int32
  acquisition_time_stamp::Int32
  physiology_time_stamp::NTuple{3,Int32}
  number_of_samples::Int16
  available_channels::Int16
  active_channels::Int16
  channel_mask::NTuple{16,Int64}
  discard_pre::Int16
  discard_post::Int16
  center_sample::Int16
  encoding_space_ref::Int16
  trajectory_dimensions::Int16
  sample_time_us::Float32
  position::NTuple{3,Float32}
  read_dir::NTuple{3,Float32}
  phase_dir::NTuple{3,Float32}
  slice_dir::NTuple{3,Float32}
  patient_table_position::NTuple{3,Float32}
  idx::EncodingCounters
  user_int::NTuple{8,Int32}
  user_float::NTuple{8,Float32}
end

function AcquisitionHeader(; version=0, flags=0, measurement_uid=0, scan_counter=0,
   acquisition_time_stamp=0, physiology_time_stamp = ntuple(i->Int32(0),3),
   number_of_samples=0, available_channels=0, active_channels=0, channel_mask=ntuple(i->Int64(0),16),
   discard_pre=0, discard_post=0, center_sample=0, encoding_space_ref=0, trajectory_dimensions=0,
   sample_time_us=0.0, position=ntuple(i->Float32(0),3), read_dir=ntuple(i->Float32(0),3),
   phase_dir=ntuple(i->Float32(0),3), slice_dir=ntuple(i->Float32(0),3),
   patient_table_position=ntuple(i->Float32(0),3), idx=EncodingCounters(),
   user_int=ntuple(i->Int32(0),8), user_float=ntuple(i->Float32(0),8))
  return AcquisitionHeader(version, flags, measurement_uid, scan_counter,
     acquisition_time_stamp, physiology_time_stamp,
     number_of_samples, available_channels, active_channels, channel_mask,
     discard_pre, discard_post, center_sample, encoding_space_ref, trajectory_dimensions,
     sample_time_us, position, read_dir, phase_dir, slice_dir, patient_table_position,
     idx, user_int, user_float)
end

mutable struct Profile
  head::AcquisitionHeader
  traj::Array{Float32,2}
  data::Array{Complex{Float32},2}
end

mutable struct RawAcquisitionData
  params::Dict{String, Any}
  profiles::Vector{Profile}
end


function trajectory(f::RawAcquisitionData)
  if f.params["trajectory"] == "cartesian"
    return trajectory("Cartesian", f.params["encodedSize"][2],
                                   f.params["encodedSize"][1])

  elseif f.params["trajectory"] == "spiral"

    encSt1 = encSteps1(f)
    encSt2 = encSteps2(f)
    sl = slices(f)
    rep = repetitions(f)

    numSl = length(unique(sl))
    numRep = length(unique(rep))
    numEncSt1 = length(unique(encSt1))
    numEncSt2 = length(unique(encSt2))

    numSampPerProfile = size(f.profiles[1].data,1)
    numChannels = size(f.profiles[1].data,2)
    D = Int(f.profiles[1].head.trajectory_dimensions)

    # remove data that should be discarded
    i1 = f.profiles[1].head.discard_pre + 1
    i2 = numSampPerProfile - f.profiles[1].head.discard_post

    traj = zeros(Float32, D, length(i1:i2), numEncSt1, numEncSt2, numSl, numRep)

    for l=1:length(f.profiles)
      traj[:, :, encSt1[l], encSt2[l], sl[l], rep[l]] .= f.profiles[l].traj[:,i1:i2]
    end

    traj_ = reshape(traj[:,:,:,:,1,1], D, :)
    return Trajectory(traj_, size(traj_,3), size(traj_,2), circular=true)
  end
  return nothing
end


function sequence(f::RawAcquisitionData)
  # TODO
end


function numChannels(f::RawAcquisitionData)
  return f.params["receiverChannels"]
end

encSteps1(f::RawAcquisitionData) =
   [f.profiles[l].head.idx.kspace_encode_step_1+1 for l=1:length(f.profiles)]
encSteps2(f::RawAcquisitionData) =
   [f.profiles[l].head.idx.kspace_encode_step_2+1 for l=1:length(f.profiles)]
slices(f::RawAcquisitionData) =
   [f.profiles[l].head.idx.slice+1 for l=1:length(f.profiles)]
repetitions(f::RawAcquisitionData) =
   [f.profiles[l].head.idx.repetition+1 for l=1:length(f.profiles)]


function rawdata(f::RawAcquisitionData)
  encSt1 = encSteps1(f)
  encSt2 = encSteps2(f)
  sl = slices(f)
  rep = repetitions(f)

  numSl = length(unique(sl))
  numRep = length(unique(rep))
  numEncSt1 = length(unique(encSt1))
  numEncSt2 = length(unique(encSt2))

  numSampPerProfile = size(f.profiles[1].data,1)
  numChan = size(f.profiles[1].data,2)

  # remove data that should be discarded
  i1 = f.profiles[1].head.discard_pre + 1
  i2 = numSampPerProfile - f.profiles[1].head.discard_post

  sorted = zeros(ComplexF32, length(i1:i2), numEncSt1, numEncSt2,
                 numChan, numSl, numRep)

  for l=1:length(f.profiles)
    sorted[:, encSt1[l], encSt2[l], :, sl[l], rep[l]] .= f.profiles[l].data[i1:i2,:]
  end

  # return map(ComplexF64, vec(sorted))
  kdata = map(ComplexF64, reshape(sorted,:,1,numChan,numRep*numSl))
  return [kdata[:,echo,:,slice] for echo=1:1,slice=1:numRep*numSl]
end

function acquisitionData(f::RawAcquisitionData)
  return AcquisitionData(trajectory(f), rawdata(f),
                          numCoils=numChannels(f),
                          numEchoes=1,
                          numSlices=length(unique(slices(f)))*length(unique(repetitions(f))),
                          encodingSize=f.params["encodedSize"] )
end
