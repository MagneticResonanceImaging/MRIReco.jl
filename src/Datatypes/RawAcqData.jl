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

function trajectory(f::RawAcquisitionData; slice::Int=1, contrast::Int=1)
  name = get(f.params, "trajectory","cartesian")
  if lowercase(name) == "cartesian"
    if f.params["encodedSize"][3]>1
      tr = trajectory("Cartesian3D", f.params["encodedSize"][2],
                                     f.params["encodedSize"][1],
                                     numSlices=f.params["encodedSize"][3])
    else
      tr = trajectory("Cartesian", f.params["encodedSize"][2],
                                   f.params["encodedSize"][1])
    end
  else
    name  = f.params["trajectory"]
    encSt1 = encSteps1(f)
    encSt2 = encSteps2(f)
    sl = slices(f)
    rep = repetitions(f)

    numSl = length(unique(sl))
    numRep = length(unique(rep))
    numEncSt1 = length(unique(encSt1))
    numEncSt2 = length(unique(encSt2))

    # assume constant number of samplings per profile
    numSampPerProfile = size(f.profiles[1].data,1)
    numChannels = size(f.profiles[1].data,2)
    D = Int(f.profiles[1].head.trajectory_dimensions)

    # remove data that should be discarded
    i1 = f.profiles[1].head.discard_pre + 1
    i2 = numSampPerProfile - f.profiles[1].head.discard_post

    traj = zeros(Float32, D, length(i1:i2), numEncSt1, numEncSt2, numSl, numRep)

    for l=1:length(f.profiles)
      if f.profiles[l].head.idx.slice+1 != slice || f.profiles[l].head.idx.contrast+1 != contrast
        continue
      end
      traj[:, :, encSt1[l], encSt2[l], sl[l], rep[l]] .= f.profiles[l].traj[:,i1:i2]
    end

    traj_ = reshape(traj[:,:,:,:,1,1], D, :)
    tr = Trajectory(traj_, size(traj_,3), size(traj_,2), circular=true)
  end
  return tr
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
contrasts(f::RawAcquisitionData) =
  [f.profiles[l].head.idx.contrast+1 for l=1:length(f.profiles)]


# get sampled points for a cartesian trajectory
function subsampleIndices(f::RawAcquisitionData; slice::Int=1, contrast::Int=1)
  idx = Int64[]
  encSt1 = encSteps1(f)
  encSt2 = encSteps2(f)
  numEncSamp, numProf, numSl = f.params["encodedSize"]
  for i=1:length(f.profiles)
    # only consider data for the specified slice, constrast and repetition
    if f.profiles[i].head.idx.slice+1 != slice || f.profiles[i].head.idx.contrast+1 != contrast || f.profiles[i].head.idx.repetition != 0
      continue
    end
    numSamp = size(f.profiles[i].data,1)  # number of measured points
    center_sample = f.profiles[i].head.center_sample+1 # idx of the center trajectory in the measurement
    tr_center_idx = div(numEncSamp,2)+1 # center of the encoded trajectory

    # start and end indices of sampled profile
    i1 = tr_center_idx - center_sample + 1 + f.profiles[i].head.discard_pre
    i2 = tr_center_idx - center_sample + numSamp - f.profiles[i].head.discard_post
    # convert to linear index
    lineIdx = collect(i1:i2) .+ numSamp*((encSt2[i]-1)*numProf + (encSt1[i]-1))
    append!(idx, lineIdx)
  end
  return sort(unique(idx))
end


function rawdata(f::RawAcquisitionData)
  encSt1 = encSteps1(f)
  encSt2 = encSteps2(f)
  sl = slices(f)
  rep = repetitions(f)
  contr = contrasts(f)

  numSl = length(unique(sl))
  numRep = length(unique(rep))
  numContr = length(unique(contr))
  numChan = size(f.profiles[1].data,2)

  kdata = Array{Array{ComplexF64,4}}(undef,numContr,numSl,numRep)
  encIdx1 = Array{Vector{Int64}}(undef,numContr,numSl)
  encIdx2 = Array{Vector{Int64}}(undef,numContr,numSl)
  numSampPerProfile = zeros(Int64, numContr,numSl)
  for j = 1:numSl
    idx_sl = findall(x->x==j,sl)
    for i = 1:numContr
      # initialize kdata with proper sized arrays for all contrasts and slices
      idx_contr = findall(x->x==i,contr)
      idx = intersect(idx_sl,idx_contr)
      numEncSt1 = length(unique(encSt1[idx]))
      numEncSt2 = length(unique(encSt2[idx]))
      numSampPerProfile[i,j] = size(f.profiles[idx[1]].data,1)
      numDiscard = f.profiles[idx[1]].head.discard_pre + f.profiles[idx[1]].head.discard_post
      for k = 1:numRep
        kdata[i,j,k] = zeros(ComplexF64,numSampPerProfile[i,j]-numDiscard,numEncSt1,numEncSt2,numChan)
      end

      # permutation from stored k-space data to the sampled locations
      encIdx1[i,j]=sort(unique(encSt1[idx]))
      encIdx2[i,j]=sort(unique(encSt2[idx]))
    end
  end

  # get k-space data
  for l=1:length(f.profiles)
    # remove data that should be discarded
    i1 = f.profiles[l].head.discard_pre + 1
    i2 = numSampPerProfile[contr[l], sl[l]] - f.profiles[l].head.discard_post
    # find k-space location idx of the current profile
    idx1 = findfirst(x->x==encSt1[l], encIdx1[contr[l],sl[l]])
    idx2 = findfirst(x->x==encSt2[l], encIdx2[contr[l],sl[l]])
    kdata[contr[l], sl[l],rep[l]][:, idx1, idx2, :] .= ComplexF64.(f.profiles[l].data[i1:i2, :])
  end

  return [reshape(kdata[c,s,r],:,numChan) for c=1:numContr, s=1:numSl, r=1:numRep]
end

function AcquisitionData(f::RawAcquisitionData)
  numSl = length(unique(slices(f)))
  numRep = length(unique(repetitions(f)))
  numContr = length(unique(contrasts(f)))
  tr = [trajectory(f,contrast=contr) for contr=1:numContr]
  subsampleIdx = [subsampleIndices(f,contrast=contr) for contr=1:numContr]
  return AcquisitionData(tr, rawdata(f),
                          numCoils=numChannels(f),
                          numEchoes=numContr,
                          numSlices=numSl,
                          numReps=numRep,
                          idx=subsampleIdx,
                          encodingSize=collect(f.params["encodedSize"]),
                          fov = collect(f.params["encodedFOV"]) )
end


function RawAcquisitionData(acqData::AcquisitionData)
  # XML header
  params = minimalHeader(Tuple(acqData.encodingSize),Tuple(acqData.fov),tr_name=string(trajectory(acqData,1)))
  # acquisition counter
  counter = 1
  # profiles
  profiles = Vector{Profile}()
  for rep = 1:acqData.numReps
    for slice=1:acqData.numSlices
      for contr = 1:acqData.numEchoes
        tr = trajectory(acqData,contr)
        profIdx = profileIdx(acqData,contr)
        numSamp = numSamplingPerProfile(tr)
        numProf = div(length(acqData.subsampleIndices[contr]),numSamp)
          for prof_tr = 1:length(profIdx)
            for slice_tr = 1:numSlices(tr)
              head = AcquisitionHeader(acqData,rep,slice,slice_tr,contr,profIdx[prof_tr],counter)
              nodes = profileNodes(tr,prof_tr,slice_tr)
              if dims(tr)==2
                data = profileData(acqData,contr,slice,rep,prof_tr)
              else
                data = profileData(acqData,contr,slice_tr,rep,prof_tr)
              end
              profile = Profile( head, Float32.(nodes), Complex{Float32}.(data) )
              push!( profiles, profile )
              counter += 1
            end
          end
      end
    end
  end

  return RawAcquisitionData(params,profiles)
end

function AcquisitionHeader(acqData::AcquisitionData, rep::Int, slice::Int, slice_tr::Int
                            , contr::Int, prof_tr::Int, counter::Int
                            ; phase::Int=1, set::Int=1, segment::Int=1, repetition::Int=1
                            , user::NTuple{8,Int16} = ntuple(i->Int16(0),8), kargs...)
  numChannels = acqData.numCoils
  tr = trajectory(acqData,contr)
  numSamples = numSamplingPerProfile(tr)
  isCartesian(tr) ? center_sample=div(numSamples,2) : center_sample=0 # needs to be fixed for partial Fourier,etc
  tr_dims = dims(tr)
  sample_time_us = acqTimePerProfile(tr)/numSamples*1.e6

  idx = EncodingCounters( Int16(prof_tr-1)
                        , Int16(slice_tr-1)
                        , Int16(rep-1)
                        , Int16(slice-1)
                        , Int16(contr-1)
                        , Int16(phase-1)
                        , Int16(repetition-1)
                        , Int16(set-1)
                        , Int16(segment-1)
                        , user )

  return AcquisitionHeader(; scan_counter=Int32(counter)
                          , number_of_samples=Int16(numSamples)
                          , available_channels=Int16(numChannels)
                          , active_channels=Int16(numChannels)
                          , center_sample=Int16(center_sample)
                          , trajectory_dimensions=Int16(tr_dims)
                          , sample_time_us=Float32(sample_time_us)
                          , read_dir=Float32.((1,0,0))
                          , phase_dir=Float32.((0,1,0))
                          , slice_dir=Float32.((0,0,1))
                          , idx=idx )
end
