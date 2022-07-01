export RawAcquisitionData, EncodingCounters, AcquisitionHeader, Profile, minimalHeader, profileIdx

"""
Encoding counters used in each Profile of a RawAcquisitionData object.
"""
Base.@kwdef mutable struct EncodingCounters
  kspace_encode_step_1::UInt16 = 0
  kspace_encode_step_2::UInt16 = 0
  average::UInt16 = 0
  slice::UInt16 = 0
  contrast::UInt16 = 0
  phase::UInt16 = 0
  repetition::UInt16 = 0
  set::UInt16 = 0
  segment::UInt16 = 0
  user::NTuple{8,UInt16} = ntuple(i->UInt16(0),8)
end

"""
Encoding counters used in each Profile of a RawAcquisitionData object.
"""
Base.@kwdef mutable struct AcquisitionHeader
  version::UInt16 = 0
  flags::UInt64 = 0
  measurement_uid::UInt32 = 0
  scan_counter::UInt32 = 0
  acquisition_time_stamp::UInt32 = 0
  physiology_time_stamp::NTuple{3,UInt32} = ntuple(i->UInt32(0),3)
  number_of_samples::UInt16 = 0
  available_channels::UInt16 = 0
  active_channels::UInt16 = 0
  channel_mask::NTuple{16,UInt64} = ntuple(i->UInt64(0),16)
  discard_pre::UInt16 = 0
  discard_post::UInt16 = 0
  center_sample::UInt16 = 0
  encoding_space_ref::UInt16 = 0
  trajectory_dimensions::UInt16 = 0
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

"""
Struct to describe data from one profile of a trajectory (from all coils).
"""
mutable struct Profile
  head::AcquisitionHeader
  traj::Array{Float32,2}
  data::Array{Complex{Float32},2}
end

"""
RawAcquisitionData object.

# Fields
* `params::Dict{String, Any}` - Dict containing the information of the XML header in ISMRMRD
* `profiles::Vector{Profile}` - Vector containing all the profiles of the acquisition

"""
mutable struct RawAcquisitionData
  params::Dict{String, Any}
  profiles::Vector{Profile}
end

"""
  trajectory(f::RawAcquisitionData; slice::Int=1, contrast::Int=1)

returns the `Trajectory` for given `slice` and `contrast` of a `RawAcquisitionData` object.
"""
function MRIBase.trajectory(f::RawAcquisitionData; slice::Int=1, contrast::Int=1)
  name = get(f.params, "trajectory","cartesian")

  sl = slices(f)
  rep = repetitions(f)
  contr = contrasts(f)

  numSl = length(unique(sl))
  numRep = length(unique(rep))
  numContr = length(unique(contr))

  encSt1 = encSteps1(f)
  encSt2 = encSteps2(f)
  numEncSt1 = length(unique(encSt1))
  numEncSt2 = length(unique(encSt2))

  # assume constant number of samplings per profile
  numSampPerProfile = size(f.profiles[1].data,1)
  numChan = size(f.profiles[1].data,2)
  D = Int(f.profiles[1].head.trajectory_dimensions)

  # remove data that should be discarded
  i1 = f.profiles[1].head.discard_pre + 1
  i2 = numSampPerProfile - f.profiles[1].head.discard_post

  if lowercase(name) == "cartesian"
    dim = ( maximum(encSteps2(f))==1 ? 2 : 3)
    if dim==3
      tr = trajectory(Float32, "Cartesian3D", f.params["encodedSize"][2],
                                     f.params["encodedSize"][1],
                                     numSlices=f.params["encodedSize"][3])
    else
      tr = trajectory(Float32, "Cartesian", f.params["encodedSize"][2],
                                   f.params["encodedSize"][1])
    end

  elseif (f.params["trajectory"]=="custom")
    numProf = Int(length(f.profiles)/(numSl * numRep * numContr))

    traj = zeros(Float32, D, length(i1:i2), numProf, 1, numSl, numRep)
    times = zeros(Float32, length(i1:i2), numProf, 1, numSl, numRep)

    counter = zeros(Int,numSl,numRep) #counter for sl and rep
    for l=1:length(f.profiles)
      if f.profiles[l].head.idx.slice+1 != slice || f.profiles[l].head.idx.contrast+1 != contrast
        continue
      end
      counter[sl[l]-minimum(sl)+1,rep[l]-minimum(rep)+1] = counter[sl[l]-minimum(sl)+1,rep[l]-minimum(rep)+1] + 1
      idx_tmp = counter[sl[l]-minimum(sl)+1,rep[l]-minimum(rep)+1]

      traj[:, :, idx_tmp, 1, sl[l]-minimum(sl)+1, rep[l]-minimum(rep)+1] .= f.profiles[l].traj[:,i1:i2]
      dt = f.profiles[l].head.sample_time_us*1e-6
      times[:, idx_tmp, 1, sl[l]-minimum(sl)+1, rep[l]-minimum(rep)+1] .= 0:dt:(length(i1:i2)-1)*dt
    end

    traj_ = reshape(traj[:,:,:,:,1,1], D, :)
    tr = Trajectory(traj_, size(traj,3), size(traj,2), circular=true, times=vec(times[:,:,:,1,1]))

  else
    name  = f.params["trajectory"]

    traj = zeros(Float32, D, length(i1:i2), numEncSt1, numEncSt2, numSl, numRep)
    times = zeros(Float32, length(i1:i2), numEncSt1, numEncSt2, numSl, numRep)

    for l=1:length(f.profiles)
      if f.profiles[l].head.idx.slice+1 != slice || f.profiles[l].head.idx.contrast+1 != contrast
        continue
      end
      traj[:, :, encSt1[l], encSt2[l], sl[l]-minimum(sl)+1, rep[l]-minimum(rep)+1] .= f.profiles[l].traj[:,i1:i2]
      dt = f.profiles[l].head.sample_time_us*1e-6
      times[:, encSt1[l], encSt2[l], sl[l]-minimum(sl)+1, rep[l]-minimum(rep)+1] .= 0:dt:(length(i1:i2)-1)*dt
    end

    traj_ = reshape(traj[:,:,:,:,1,1], D, :)
    # tr = Trajectory(traj_, size(traj_,3), size(traj_,2), circular=true)
    tr = Trajectory(traj_, size(traj,3), size(traj,2), circular=true, times=vec(times[:,:,:,1,1]))
  end
  return tr
end

function sequence(f::RawAcquisitionData)
  # TODO
end

"""
    numChannels(f::RawAcquisitionData)

returns the number of channels in a `RawAcquisitionData` object.
"""
function numChannels(f::RawAcquisitionData)
  return f.profiles[1].head.active_channels
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

  """
    subsampleIndices(f::RawAcquisitionData; slice::Int=1, contrast::Int=1)

  returns the sampled indices for a given `slice` and `contrast` in a `RawAcquisitionData` object.
  """
function subsampleIndices(f::RawAcquisitionData; slice::Int=1, contrast::Int=1, estimateProfileCenter::Bool=false)
  idx = Int64[]
  encSt1 = encSteps1(f)
  encSt2 = encSteps2(f)
  numEncSamp, numProf, numSl = f.params["encodedSize"]
  for i=1:length(f.profiles)
    # only consider data for the specified slice, contrast and repetition
    if f.profiles[i].head.idx.slice+1 != slice || f.profiles[i].head.idx.contrast+1 != contrast || f.profiles[i].head.idx.repetition != 0
      continue
    end
    numSamp = size(f.profiles[i].data,1)  # number of measured points
    if estimateProfileCenter
      center_sample = div(numSamp,2)+1
    else
      center_sample = f.profiles[i].head.center_sample+1 # idx of the center trajectory in the measurement
    end
    tr_center_idx = div(numEncSamp,2)+1 # center of the encoded trajectory

    # start and end indices of sampled profile
    i1 = tr_center_idx - center_sample + 1 + f.profiles[i].head.discard_pre
    i2 = tr_center_idx - center_sample + numSamp - f.profiles[i].head.discard_post
    # convert to linear index
    lineIdx = collect(i1:i2) .+ numSamp*((encSt2[i]-1)*numProf + (encSt1[i]-1))
    append!(idx, lineIdx)
  end

  return unique(idx)
end

"""
    rawdata(f::RawAcquisitionData)

returns the rawdata contained `RawAcquisitionData` object.
The output is an `Array{Matrix{ComplexF64},3}`, which can be stored in a `AcquisitionData` object.
"""
function rawdata(f::RawAcquisitionData; slice::Int=1, contrast::Int=1, repetition::Int=1)

  # find profiles corresponding to the given slice contrast and repitiion
  sl = slices(f)
  rep = repetitions(f)
  contr = contrasts(f)
  idx_sl = findall(x->x==slice,sl)
  idx_contr = findall(x->x==contrast,contr)
  idx_rep = findall(x->x==repetition,rep)
  idx = intersect(idx_sl,idx_contr,idx_rep)

  numSl = length(unique(sl))
  numRep = length(unique(rep))
  numContr = length(unique(contr))

  # number of unique combination of encoding statuses
  encSt1 = encSteps1(f)[idx]
  encSt2 = encSteps2(f)[idx]
  # find first occurence of each combination of encoding statuses (as `unique` returns)
  idx_unique = uniqueidx(hcat(encSt1,encSt2))
  numEncSt = length(idx_unique)

  # allocate space for k-space data
  # assume the same number of samples for all profiles
  numSampPerProfile, numChan = size(f.profiles[idx[1]].data)
  numSampPerProfile -= (f.profiles[idx[1]].head.discard_pre+f.profiles[idx[1]].head.discard_post)
  
  if f.params["trajectory"] == "custom"
    numProf = Int(length(f.profiles)/(numSl * numRep * numContr))
    kdata = zeros(typeof(f.profiles[1].data[1, 1]), numSampPerProfile, numProf, numChan)
    posIdx = 1:length(f.profiles)
  else # cartesian case
    kdata = zeros(typeof(f.profiles[1].data[1, 1]), numSampPerProfile, numEncSt, numChan)
    posIdx = idx[idx_unique]
  end

  # store one profile (kspace data) for each unique encoding status
  cnt = 1
  for l in posIdx # TODO: set of profiles (with unique encoding status)
    if f.profiles[l].head.idx.slice+1 != slice || f.profiles[l].head.idx.contrast+1 != contrast || f.profiles[l].head.idx.repetition+1 != repetition
      continue
    end
    i1 = f.profiles[l].head.discard_pre + 1
    i2 = i1+numSampPerProfile-1
    kdata[:,cnt,:] .= f.profiles[l].data[i1:i2, :]
    cnt += 1
  end
  
  return reshape(kdata, :, numChan)
end


"""
    AcquisitionData(f::RawAcquisitionData; estimateProfileCenter::Bool=false)

converts `RawAcquisitionData` into the equivalent `AcquisitionData` object.
"""
function AcquisitionData(f::RawAcquisitionData; estimateProfileCenter::Bool=false, OffsetBruker=false)
  contrs = sort(unique(contrasts(f)))
  sls = sort(unique(slices(f)))
  reps = sort(unique(repetitions(f)))
  numContr = length(unique(contrasts(f)))
  numSl = length(unique(slices(f)))
  numRep = length(unique(repetitions(f)))
  # tr = [trajectory(f,contrast=contr) for contr=1:numContr]
  # subsampleIdx = [subsampleIndices(f,contrast=contr,estimateProfileCenter=estimateProfileCenter) for contr=1:numContr]
  # kdata = [rawdata(f; contrast=i, slice=j, repetition=k) for i=1:numContr, j=1:numSl, k=1:numRep]

  
  tr = [trajectory(f,contrast=contr) for contr=contrs]
  subsampleIdx = [subsampleIndices(f,contrast=contr,estimateProfileCenter=estimateProfileCenter) for contr=contrs]
  kdata = [rawdata(f; contrast=i, slice=j, repetition=k) for i=contrs, j=sls, k=reps]

  if OffsetBruker
    @info "Offset correction is performed"
    for i = 1:length(tr) #vector of trajectory for contrasts
      for k in reps
        ROT = [[f.profiles[1].head.read_dir...] [f.profiles[1].head.phase_dir...] [f.profiles[1].head.slice_dir...]]
        shift = inv(ROT) * [f.profiles[1].head.position...]./f.params["encodedFOV"]
        shift = shift .*[0 1 -1]'; @info "read offset not performed (only phase + slice)";
        shift = shift .* f.params["encodedSize"]/2
        phase_nodes = exp.(-2π * im * tr[i].nodes' * shift)
        
        [kdata[i,j,k] = kdata[i,j,k] .* phase_nodes for i=contrs, j=sls]
      end
    end
  end

  return AcquisitionData(tr, kdata,
                          idx=subsampleIdx,
                          encodingSize=collect(f.params["encodedSize"]),
                          fov = collect(f.params["encodedFOV"]) )
end

"""
    RawAcquisitionData(acqData::AcquisitionData)

converts `acqData` into the equivalent `RawAcquisitionData` object.
"""
function RawAcquisitionData(acqData::AcquisitionData)
  # XML header
  params = minimalHeader(Tuple(acqData.encodingSize),Tuple(acqData.fov),tr_name=string(trajectory(acqData,1)))
  # acquisition counter
  counter = 1
  # profiles
  profiles = Vector{Profile}()
  for rep = 1:numRepetitions(acqData)
    for slice=1:numSlices(acqData)
      for contr = 1:numContrasts(acqData)
        tr = trajectory(acqData,contr)
        profIdx = profileIdx(acqData,contr)
        numSamp = numSamplingPerProfile(tr)
        numProf = div(length(acqData.subsampleIndices[contr]),numSamp)
          for prof_tr = 1:length(profIdx)
            for slice_tr = 1:numSlices(tr)
              head = AcquisitionHeader(acqData,rep,slice,slice_tr,contr,profIdx[prof_tr],counter)
              nodes = profileNodes(tr,prof_tr,slice_tr)
              if ndims(tr)==2
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
  numChan = numChannels(acqData)
  tr = trajectory(acqData,contr)
  numSamples = numSamplingPerProfile(tr)
  isCartesian(tr) ? center_sample=div(numSamples,2) : center_sample=0 # needs to be fixed for partial Fourier,etc
  tr_dims = ndims(tr)
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
                          , available_channels=Int16(numChan)
                          , active_channels=Int16(numChan)
                          , center_sample=Int16(center_sample)
                          , trajectory_dimensions=Int16(tr_dims)
                          , sample_time_us=Float32(sample_time_us)
                          , read_dir=Float32.((1,0,0))
                          , phase_dir=Float32.((0,1,0))
                          , slice_dir=Float32.((0,0,1))
                          , idx=idx )
end

# get indices of the first occurence of unique elements in the matrix x
function uniqueidx(x::Matrix{T}) where T
  uniqueset = Set{Vector{T}}()
  ex = eachindex(x[:,1])
  idxs = Vector{eltype(ex)}()
  for i in ex
      xi = x[i,:]
      if !(xi in uniqueset)
          push!(idxs, i)
          push!(uniqueset, xi)
      end
  end
  idxs
end

function minimalHeader(encodingSize::NTuple{3,Int},fov::NTuple{3,AbstractFloat};f_res::Integer=1,tr_name::AbstractString="cartesian",numChannels::Int=1)
  params = Dict{String,Any}()
  params["H1resonanceFrequency_Hz"] = f_res
  params["encodedSize"] = collect(encodingSize)
  params["encodedFOV"] = collect(fov)
  params["trajectory"] = tr_name
  params["receiverChannels"] = numChannels

  return params
end

function profileIdx(acqData::AcquisitionData,contr::Int)
  tr = trajectory(acqData,contr)
  numSamp = numSamplingPerProfile(tr)
  numProf = div(length(acqData.subsampleIndices[contr]),numSamp)
  idx = [div(acqData.subsampleIndices[contr][numSamp*(prof-1)+1],numSamp)+1 for prof=1:numProf]
end