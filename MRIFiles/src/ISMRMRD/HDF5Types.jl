import HDF5.API: h5t_create, h5t_insert, h5t_array_create, h5t_close,
             h5t_vlen_create, hsize_t, h5s_create_simple, h5d_create, h5d_write,
             h5d_close, h5s_close, h5p_create, h5p_set_chunk, h5t_get_member_offset,
             h5t_get_member_type, h5t_get_member_offset, h5t_get_member_type,
             H5T_COMPOUND, H5P_DEFAULT, H5P_DATASET_CREATE, H5S_ALL

import HDF5: hdf5_type_id

const ISMRMRD_HEADER_SIZE = 340
const ISMRMRD_HEADER_PLUS_POINTERS_SIZE = 372
const ENCODING_COUNTERS_SIZE = 34

function get_hdf5type_complex(::Type{T}=Float32) where T
  datatype = h5t_create(H5T_COMPOUND,2*sizeof(T))
  h5t_insert(datatype, "real", 0 , hdf5_type_id(T))
  h5t_insert(datatype, "imag", sizeof(T) , hdf5_type_id(T))
  #HDF5.h5t_close(datatype)
  return datatype
end

#=
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
=#

function get_hdf5type_encoding()
  datatype = h5t_create(H5T_COMPOUND, ENCODING_COUNTERS_SIZE )
  h5t_insert(datatype, "kspace_encode_step_1", 0, hdf5_type_id(UInt16))
  h5t_insert(datatype, "kspace_encode_step_2", 2 , hdf5_type_id(UInt16))
  h5t_insert(datatype, "average", 4, hdf5_type_id(UInt16))
  h5t_insert(datatype, "slice", 6, hdf5_type_id(UInt16))
  h5t_insert(datatype, "contrast", 8, hdf5_type_id(UInt16))
  h5t_insert(datatype, "phase", 10, hdf5_type_id(UInt16))
  h5t_insert(datatype, "repetition", 12, hdf5_type_id(UInt16))
  h5t_insert(datatype, "set", 14, hdf5_type_id(UInt16))
  h5t_insert(datatype, "segment", 16, hdf5_type_id(UInt16))

  hsz = hsize_t[8]
  d = h5t_array_create(hdf5_type_id(UInt16), convert(Cuint, 1), hsz)
  h5t_insert(datatype, "user", 18, d)
  h5t_close(d)
  return datatype
end

#=

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
=#


function get_hdf5type_acquisitionheader()
  datatype = h5t_create(H5T_COMPOUND, ISMRMRD_HEADER_SIZE )
  h5t_insert(datatype, "version", 0, hdf5_type_id(UInt16))
  h5t_insert(datatype, "flags", 2 , hdf5_type_id(UInt64))
  h5t_insert(datatype, "measurement_uid", 10, hdf5_type_id(UInt32))
  h5t_insert(datatype, "scan_counter", 14, hdf5_type_id(UInt32))
  h5t_insert(datatype, "acquisition_time_stamp", 18, hdf5_type_id(UInt32))

  hsz1 = hsize_t[3]
  d1 = h5t_array_create(hdf5_type_id(UInt32), convert(Cuint, 1), hsz1)
  h5t_insert(datatype, "physiology_time_stamp", 22, d1)
  h5t_close(d1)

  h5t_insert(datatype, "number_of_samples", 34, hdf5_type_id(UInt16))
  h5t_insert(datatype, "available_channels", 36, hdf5_type_id(UInt16))
  h5t_insert(datatype, "active_channels", 38, hdf5_type_id(UInt16))

  hsz2 = hsize_t[16]
  d2 = h5t_array_create(hdf5_type_id(UInt64), convert(Cuint, 1), hsz2)
  h5t_insert(datatype, "channel_mask", 40, d2)
  h5t_close(d2)

  h5t_insert(datatype, "discard_pre", 168, hdf5_type_id(UInt16))
  h5t_insert(datatype, "discard_post", 170, hdf5_type_id(UInt16))
  h5t_insert(datatype, "center_sample", 172, hdf5_type_id(UInt16))
  h5t_insert(datatype, "encoding_space_ref", 174, hdf5_type_id(UInt16))
  h5t_insert(datatype, "trajectory_dimensions", 176, hdf5_type_id(UInt16))
  h5t_insert(datatype, "sample_time_us", 178, hdf5_type_id(Float32))

  hsz3 = hsize_t[3]
  d3 = h5t_array_create(hdf5_type_id(Float32), convert(Cuint, 1), hsz3)
  h5t_insert(datatype, "position", 182, d3)
  h5t_insert(datatype, "read_dir", 194, d3)
  h5t_insert(datatype, "phase_dir", 206, d3)
  h5t_insert(datatype, "slice_dir", 218, d3)
  h5t_insert(datatype, "patient_table_position", 230, d3)
  h5t_close(d3)

  denc = get_hdf5type_encoding()
  h5t_insert(datatype, "idx", 242, denc)
  h5t_close(denc)

  hsz4 = hsize_t[8]
  d4 = h5t_array_create(hdf5_type_id(Int32), convert(Cuint, 1), hsz4)
  h5t_insert(datatype, "user_int", 276, d4)
  h5t_close(d4)

  hsz5 = hsize_t[8]
  d5 = h5t_array_create(hdf5_type_id(Float32), convert(Cuint, 1), hsz5)
  h5t_insert(datatype, "user_float", 308, d5)
  h5t_close(d5)

  #HDF5.h5t_close(datatype)
  return datatype
end

#=
typedef struct HDF5_Acquisition
{
    ISMRMRD_AcquisitionHeader head;
    hvl_t traj;
    hvl_t data;
} HDF5_Acquisition;
=#

struct HVL_T{T}
  len::Int64
  p::Ptr{T}
end

function get_hdf5type_acquisition()
  datatype = h5t_create(H5T_COMPOUND, ISMRMRD_HEADER_PLUS_POINTERS_SIZE )
  d1 = get_hdf5type_acquisitionheader()
  h5t_insert(datatype, "head", 0, d1)
  h5t_close(d1)

  vlvartype = h5t_vlen_create(hdf5_type_id(Float32))
  h5t_insert(datatype, "traj", 340, vlvartype)
  h5t_close(vlvartype)

  #vartype = get_hdf5type_complex(Float32)
  vartype = hdf5_type_id(Float32) #for some reason this is stored as float
  vlvartype = h5t_vlen_create(vartype)
  h5t_insert(datatype, "data", 356, vlvartype)
  h5t_close(vlvartype)
  #h5t_close(vartype)

  return datatype
end

function get_hdf5type_acquisition_header_only()
  datatype = h5t_create(H5T_COMPOUND, ISMRMRD_HEADER_SIZE )
  d1 = get_hdf5type_acquisitionheader()
  h5t_insert(datatype, "head", 0, d1)
  h5t_close(d1)

  return datatype
end

# Next we define several functions for converting the structs into (packed)
# binary blobs. We use IOBuffers for that which internally encode
# the data as UInt8 vectors


function Base.write(io::IOBuffer, tup::NTuple{N, T}) where {N,T}
  for d in tup
    write(io, d)
  end
end

function Base.write(io::IOBuffer, acq::EncodingCounters)
  for key ∈ propertynames(acq)
    write(io, getfield(acq, key))
  end
end

function Base.write(io::IOBuffer, acq::AcquisitionHeader)
  for key ∈ propertynames(acq)
    write(io, getfield(acq, key))
  end
end

function Base.read(io::IOBuffer, ::Type{NTuple{D, T}}) where {D,T}
  return ntuple(d->read(io, T), D)
end
    
function Base.read(io::IOBuffer, ::Type{EncodingCounters})
  enc = EncodingCounters()
  for key ∈ propertynames(enc)
    d = read(io, fieldtype(EncodingCounters, key))
    setfield!(enc, key, d)
  end
  return enc
end

function Base.read(io::IOBuffer, ::Type{AcquisitionHeader})
  acq = AcquisitionHeader()
  for key ∈ propertynames(acq)
    d = read(io, fieldtype(AcquisitionHeader, key))
    setfield!(acq, key, d)
  end
  return acq
end

function Base.write(io::IOBuffer, d::HVL_T{T}) where T
  write(io, d.len)
  write(io, d.p)
end




function writeProfiles(file, dataset, profiles::Vector{Profile})

  profiles_hdf5 = Vector{NTuple{ISMRMRD_HEADER_PLUS_POINTERS_SIZE, UInt8}}(undef,length(profiles))

  for l=1:length(profiles)

    t = HVL_T(length(profiles[l].traj), pointer(profiles[l].traj))
    #for some reason this is stored as float
    d = HVL_T(2*length(profiles[l].data), pointer(reinterpret(Float32,profiles[l].data)))

    iobuf = IOBuffer()
    write(iobuf, profiles[l].head )
    write(iobuf, t )
    write(iobuf, d)
    profiles_hdf5[l] = Tuple(iobuf.data)
  end

  d_type_compound = get_hdf5type_acquisition()

  shape = hsize_t[length(profiles)]
  maxshape = hsize_t[typemax(hsize_t)]
  chunkshape = hsize_t[1]
  space = h5s_create_simple(1, shape, maxshape)

  props = h5p_create(H5P_DATASET_CREATE)
  # enable chunking so that the dataset is extensible
  h5status = h5p_set_chunk(props, 1, chunkshape)

  dset_compound = h5d_create(file, dataset, d_type_compound, space,
                                  H5P_DEFAULT, props, H5P_DEFAULT)
  h5s_close(space)
  h5d_write(dset_compound, d_type_compound, H5S_ALL, H5S_ALL, H5P_DEFAULT, profiles_hdf5)
  h5d_close(dset_compound)
  h5t_close(d_type_compound)
end
