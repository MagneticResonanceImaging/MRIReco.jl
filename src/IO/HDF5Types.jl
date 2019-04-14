import HDF5: h5t_create, hdf5_type_id, h5t_insert, h5t_array_create, h5t_close,
             h5t_vlen_create, Hsize, h5s_create_simple, h5d_create, h5d_write,
             h5d_close, h5s_close, h5p_create, h5p_set_chunk, h5t_get_member_offset,
             h5t_get_member_type, h5t_get_member_offset, h5t_get_member_type



function get_hdf5type_complex(::Type{T}=Float32) where T
  datatype = HDF5.h5t_create(HDF5.H5T_COMPOUND,2*sizeof(T))
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
  off = i -> Cint(fieldoffset(EncodingCounters,i))
  datatype = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(EncodingCounters) )
  h5t_insert(datatype, "kspace_encode_step_1", off(1) , hdf5_type_id(UInt16))
  h5t_insert(datatype, "kspace_encode_step_2", off(2) , hdf5_type_id(UInt16))
  h5t_insert(datatype, "average", off(3), hdf5_type_id(UInt16))
  h5t_insert(datatype, "slice", off(4), hdf5_type_id(UInt16))
  h5t_insert(datatype, "contrast", off(5), hdf5_type_id(UInt16))
  h5t_insert(datatype, "phase", off(6), hdf5_type_id(UInt16))
  h5t_insert(datatype, "repetition", off(7), hdf5_type_id(UInt16))
  h5t_insert(datatype, "set", off(8), hdf5_type_id(UInt16))
  h5t_insert(datatype, "segment", off(9), hdf5_type_id(UInt16))

  hsz = Hsize[8]
  d = h5t_array_create(hdf5_type_id(UInt16), convert(Cuint, 1), hsz)
  h5t_insert(datatype, "user", off(10), d)
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
  off = i -> Cint(fieldoffset(AcquisitionHeader,i))
  datatype = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(AcquisitionHeader) )
  h5t_insert(datatype, "version", off(1) , hdf5_type_id(UInt16))
  h5t_insert(datatype, "flags", off(2) , hdf5_type_id(UInt64))
  h5t_insert(datatype, "measurement_uid", off(3), hdf5_type_id(UInt32))
  h5t_insert(datatype, "scan_counter", off(4), hdf5_type_id(UInt32))
  h5t_insert(datatype, "acquisition_time_stamp", off(5), hdf5_type_id(UInt32))

  hsz1 = Hsize[3]
  d1 = h5t_array_create(hdf5_type_id(UInt32), convert(Cuint, 1), hsz1)
  h5t_insert(datatype, "physiology_time_stamp", off(6), d1)
  h5t_close(d1)

  h5t_insert(datatype, "number_of_samples", off(7), hdf5_type_id(UInt16))
  h5t_insert(datatype, "available_channels", off(8), hdf5_type_id(UInt16))
  h5t_insert(datatype, "active_channels", off(9), hdf5_type_id(UInt16))

  hsz2 = Hsize[16]
  d2 = h5t_array_create(hdf5_type_id(UInt64), convert(Cuint, 1), hsz2)
  h5t_insert(datatype, "channel_mask", off(10), d2)
  h5t_close(d2)

  h5t_insert(datatype, "discard_pre", off(11), hdf5_type_id(UInt16))
  h5t_insert(datatype, "discard_post", off(12), hdf5_type_id(UInt16))
  h5t_insert(datatype, "center_sample", off(13), hdf5_type_id(UInt16))
  h5t_insert(datatype, "encoding_space_ref", off(14), hdf5_type_id(UInt16))
  h5t_insert(datatype, "trajectory_dimensions", off(15), hdf5_type_id(UInt16))
  h5t_insert(datatype, "sample_time_us", off(16), hdf5_type_id(Float32))

  hsz3 = Hsize[3]
  d3 = h5t_array_create(hdf5_type_id(Float32), convert(Cuint, 1), hsz3)
  h5t_insert(datatype, "position", off(17), d3)
  h5t_insert(datatype, "read_dir", off(18), d3)
  h5t_insert(datatype, "phase_dir", off(19), d3)
  h5t_insert(datatype, "slice_dir", off(20), d3)
  h5t_insert(datatype, "patient_table_position", off(21), d3)
  h5t_close(d3)

  denc = get_hdf5type_encoding()
  h5t_insert(datatype, "idx", off(22), denc)
  h5t_close(denc)

  hsz4 = Hsize[8]
  d4 = h5t_array_create(hdf5_type_id(Int32), convert(Cuint, 1), hsz4)
  h5t_insert(datatype, "user_int", off(23), d4)
  h5t_close(d4)

  hsz5 = Hsize[8]
  d5 = h5t_array_create(hdf5_type_id(Float32), convert(Cuint, 1), hsz5)
  h5t_insert(datatype, "user_float", off(24), d5)
  h5t_close(d5)

  #HDF5.h5t_close(datatype)
  return datatype
end

#=
typedef struct HDF5_Acquisiton
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

struct HDF5_Acquisition
  head::AcquisitionHeader
  traj::HVL_T{Float32}
  data::HVL_T{Float32} #for some reason this is stored as float
end

function get_hdf5type_acquisition()
  off = i -> Cint(fieldoffset(HDF5_Acquisition,i))
  datatype = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(HDF5_Acquisition) )
  d1 = get_hdf5type_acquisitionheader()
  h5t_insert(datatype, "head", off(1), d1)
  h5t_close(d1)

  vlvartype = h5t_vlen_create(hdf5_type_id(Float32))
  h5t_insert(datatype, "traj", off(2), vlvartype)
  h5t_close(vlvartype)

  #vartype = get_hdf5type_complex(Float32)
  vartype = hdf5_type_id(Float32) #for some reason this is stored as float
  vlvartype = h5t_vlen_create(vartype)
  h5t_insert(datatype, "data", off(3), vlvartype)
  h5t_close(vlvartype)
  #h5t_close(vartype)

  return datatype
end


function writeProfiles(file, dataset, profiles::Vector{Profile})

  profiles_hdf5 = Vector{HDF5_Acquisition}(undef,length(profiles))
  for l=1:length(profiles)

    t = HVL_T(length(profiles[l].traj), pointer(profiles[l].traj))
    #for some reason this is stored as float
    d = HVL_T(2*length(profiles[l].data), pointer(reinterpret(Float32,profiles[l].data)))

    profiles_hdf5[l] = HDF5_Acquisition(profiles[l].head,t,d)
  end


  d_type_compound = get_hdf5type_acquisition()

  shape = Hsize[length(profiles)]
  maxshape = Hsize[typemax(Hsize)]
  chunkshape = Hsize[1]#length(profiles)]
  space = h5s_create_simple(1, shape, maxshape)

  props = h5p_create(HDF5.H5P_DATASET_CREATE)
  # enable chunking so that the dataset is extensible
  h5status = h5p_set_chunk(props, 1, chunkshape)

  dset_compound = h5d_create(file, dataset, d_type_compound, space,
                                  HDF5.H5P_DEFAULT, props, HDF5.H5P_DEFAULT)
  h5s_close(space)
  h5d_write(dset_compound, d_type_compound, HDF5.H5S_ALL, HDF5.H5S_ALL, HDF5.H5P_DEFAULT, profiles_hdf5)
  h5d_close(dset_compound)
  h5t_close(d_type_compound)
end
