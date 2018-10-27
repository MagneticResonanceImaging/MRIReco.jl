export ISMRMRD, ISMRMRDEncodingCounters, ISMRMRDAcquisitionHeader, numRepetitions

struct ISMRMRDEncodingCounters
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

struct ISMRMRDAcquisitionHeader
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
  idx::ISMRMRDEncodingCounters
  user_int::Vector{Int32}
  user_float::Vector{Float32}
end


struct ISMRMRD <: MRIFile
  filename::String
  params::Dict
  head::Array{ISMRMRDAcquisitionHeader,1}
  traj::Array{Float32,2}
  data::Array{Complex{Float32},3}
  xdoc
end

function ISMRMRD(filename::String)
  headerStr = h5read(filename, "/dataset/xml")
  xdoc = parse_string(headerStr[1])

  header = parse_xml_header(xdoc)

  d = h5read(filename, "/dataset/data")

  M = length(d)

  head = ISMRMRDAcquisitionHeader[]
  traj = Float32[]
  data = Any[]

  chan = header["receiverChannels"]

  for m=1:M

    push!(head, read_header(d[m].data[1]))

    N = reinterpret(Int64, d[m].data[2][1:8])[1]
    if N > 0
      ptr = reinterpret(Int64, d[m].data[2][9:end])[1]
      U = unsafe_wrap(Array,Ptr{Float32}(ptr),(N,))
      append!(traj, U)
    end

    N = reinterpret(Int64, d[m].data[3][1:8])[1]
    if N > 0
      ptr = reinterpret(Int64, d[m].data[3][9:end])[1]
      U = unsafe_wrap(Array,Ptr{ComplexF32}(ptr),(div(N,2),))
      push!(data, reshape(U,div(N,2*chan),chan))
    end
  end

  dataNew = zeros(ComplexF32, size(data[1],1), size(data[1],2), length(data) )
  for m=1:M
    dataNew[:,:,m] = data[m]
  end

  if !isempty(traj)
    trajNew = reshape(traj,:,M)
  else
    trajNew = zeros(Float32,0,0)
  end

  return ISMRMRD(filename, header, head, trajNew, dataNew, xdoc)
end

function parse_xml_header(xdoc)
  header = Dict{String,Any}()

  e = get_elements_by_tagname(LightXML.root(xdoc),"acquisitionSystemInformation")[1]

  header["receiverChannels"] = parse(Int,content(get_elements_by_tagname(e,"receiverChannels")[1]))

  e = get_elements_by_tagname(LightXML.root(xdoc),"encoding")[1]
  header["trajectory"] = content(get_elements_by_tagname(e,"trajectory")[1])

  e = LightXML.root(xdoc)["encoding"][1]["encodedSpace"][1]["matrixSize"][1]
  header["encodedMatrixSize"] = [parse(Int,content(e["x"][1])),
                          parse(Int,content(e["y"][1])),
                          parse(Int,content(e["z"][1]))]

  return header
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
  idx = ISMRMRDEncodingCounters(kspace_encode_step_1,kspace_encode_step_2,
                      average,slice,contrast,phase,repetition,set,segment,user)

  user_int = read!(buf,zeros(Int32,8))
  user_float = read!(buf,zeros(Float32,8))

  head = ISMRMRDAcquisitionHeader( version, flags, measurement_uid,
    scan_counter, acquisition_time_stamp, physiology_time_stamp, number_of_samples,
    available_channels, active_channels, channel_mask, discard_pre, discard_post,
    center_sample, encoding_space_ref, trajectory_dimensions,sample_time_us,
    position, read_dir, phase_dir, slice_dir, patient_table_position, idx, user_int,
    user_float
    )

  return head
end


function trajectory(f::ISMRMRD)
  if f.params["trajectory"] == "cartesian"
    return trajectory("Cartesian", f.params["encodedMatrixSize"][2],
                                   f.params["encodedMatrixSize"][1])
  elseif f.params["trajectory"] == "spiral"


    #=function trajectory(f::DFFile)
      numSamplingPerProfile, numProfiles, nodes = open(f.trajfilename,"r") do fd
        tmp1,numSamplingPerProfile,numProfiles,tmp2= read!(fd,Array{Int32}(undef,4))
        nodes = read!(fd, Array{Float32}(undef, 2, numSamplingPerProfile,numProfiles))
        return numSamplingPerProfile, numProfiles, nodes
      end

      return CustomTrajectory(numProfiles, numSamplingPerProfile, vec(nodes))
    end=#


    return nothing
  end
  return nothing
end


function sequence(f::ISMRMRD)
  # TODO
end

function numRepetitions(f::ISMRMRD)
  f.head[end].idx.repetition - f.head[1].idx.repetition + 1
end

function numChannels(f::ISMRMRD)
  return f.params["receiverChannels"]
end

function findIndices(f::ISMRMRD, repetition=1, slice=1)
  idx = zeros(Int,f.params["encodedMatrixSize"][2])
  i = 1
  for l=1:length(f.head)
    if f.head[l].idx.slice+1 == slice &&
       f.head[l].idx.repetition+1 == repetition
      idx[ f.head[l].idx.kspace_encode_step_1+1 ] = l
    end
  end
  return idx
end

function rawdata(f::ISMRMRD, repetition=1)
  return map(ComplexF64,vec(permutedims(f.data[:,:,findIndices(f,repetition)],(1,3,2))))
end

function acquisitionData(f::ISMRMRD)
  return AcquisitionData(trajectory(f), rawdata(f),
                          numCoils=numChannels(f), numEchoes=1, numSlices=1)
end
