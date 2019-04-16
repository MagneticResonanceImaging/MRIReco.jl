export RawAcquisitionData, EncodingCounters, AcquisitionHeader

Base.@kwdef struct EncodingCounters
  kspace_encode_step_1::Int16 = 0
  kspace_encode_step_2::Int16 = 0
  average::Int16 = 0
  slice::Int16 = 0
  contrast::Int16 = 0
  phase::Int16 = 0
  repetition::Int16 = 0
  set::Int16 = 0
  segment::Int16 = 0
  user::NTuple{8,Int16} = ntuple(i->Int16(0),8)
end

Base.@kwdef struct AcquisitionHeader
  version::Int16 = 0
  flags::Int64 = 0
  measurement_uid::Int32 = 0
  scan_counter::Int32 = 0
  acquisition_time_stamp::Int32 = 0
  physiology_time_stamp::NTuple{3,Int32} = ntuple(i->Int32(0),3)
  number_of_samples::Int16 = 0
  available_channels::Int16 = 0
  active_channels::Int16 = 0
  channel_mask::NTuple{16,Int64} = ntuple(i->Int64(0),16)
  discard_pre::Int16 = 0
  discard_post::Int16 = 0
  center_sample::Int16 = 0
  encoding_space_ref::Int16 = 0
  trajectory_dimensions::Int16 = 0
  sample_time_us::Float32 = 0.0
  position::NTuple{3,Float32} = ntuple(i->Float32(0),3)
  read_dir::NTuple{3,Float32} = ntuple(i->Float32(0),3)
  phase_dir::NTuple{3,Float32} = ntuple(i->Float32(0),3)
  slice_dir::NTuple{3,Float32} = ntuple(i->Float32(0),3)
  patient_table_position::NTuple{3,Float32} = ntuple(i->Float32(0),3)
  idx::EncodingCounters = EncodingCounters()
  user_int::NTuple{8,Int32} = ntuple(i->Int32(0),8)
  user_float::NTuple{8,Float32} = ntuple(i->Float32(0),8)
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
