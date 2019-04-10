export saveasIBIFile, trajectory, sequence, rawdata, acquisitionData,
       loadIBIFile, loadParams, saveParams

using HDF5

export MRIFileIBI

# Only hold the handle to the file an load/save all data on demand.
struct MRIFileIBI <: MRIFile
  filename::String
end

### reading ###

function trajectory(f::MRIFileIBI)
  fid = h5open(f.filename, "r")
  tr = fid["trajectory"]

  params = Dict{Symbol, Any}()
  for field = names(tr)
    params[Symbol(field)] = read(tr[field])
  end

  close(fid)
  return trajectory(params[:name], params[:numProfiles], params[:numSamplingPerProfile]; params...)
end

function sequenceInfo(f::MRIFileIBI)
  fid = h5open(f.filename, "r")

  # general sequence inforamation
  seq = fid["sequence"]
  seqParams = Dict{Symbol, Any}()
  for field = names(seq)
    seqParams[Symbol(field)] = read(seq[field])
  end
  close(fid)

  return seqParams
end

function rawdata(f::MRIFileIBI)
  data = h5read(f.filename, "/rawdata")
  # workaround for hdf5 not supporting complex
  return collect(reshape(reinterpret(Complex{eltype(data)}, vec(data)), (size(data)[2:end]...,) ))
end

function acquisitionData(f::MRIFileIBI)
  return AcquisitionData(sequenceInfo(f),
                        rawdata(f),
                        h5read(f.filename, "numEchoes"),
                        h5read(f.filename, "numCoils"),
                        h5read(f.filename, "numSlices"),
                        h5read(f.filename, "samplePointer"),
                        h5read(f.filename, "samplingIdx"),
                        h5read(f.filename, "encodingSize"),
                        h5read(f.filename, "fov"))
end

### writing ###

function saveasIBIFile(filename::AbstractString, rawdata::Array{Complex{T}}, tr::Trajectory) where T<:Real


  h5open(filename, "w") do file

    # Trajectory Information
    write(file, "trajectory/name", string(tr))
    write(file, "trajectory/numProfiles", tr.numProfiles)
    write(file, "trajectory/numSamplingPerProfile", tr.numSamplingPerProfile)
    write(file, "trajectory/TE", tr.TE)
    write(file, "trajectory/AQ", tr.AQ)

    rawdata_real = collect(reshape(reinterpret(T, vec(rawdata)), (2,size(rawdata)...)))
    write(file, "/rawdata", rawdata_real)

  end
end

function saveasIBIFile(filename::AbstractString, acqData::AcquisitionData)

  h5open(filename, "w") do file

    # Sequence Information
    seqInfo = sequenceInfo(seq)
    for key in keys(seqInfo)
      a = get(seqInfo,key,"")
      write(file, "sequence/"*string(k), a)
    end

    write(file, "trajectory/name", string(trajectory(acqData)))
    tr = trajectory(acqData)
    for field in fieldnames(typeof(tr))
      a = getfield(tr,field)
      write( file, "trajectory/"*string(field), a )
    end

    write(file, "numEchoes", numEchoes(acqData))

    # number of Coils, and sampled points and number of slices
    write(file, "numCoils", acqData.numCoils)
    write(file, "samplingIdx", acqData.subsampleIndices)
    write(file, "numSlices", acqData.numSlices)

    # pointer to the data corresponding to a given echo, coil and slice
    write(file, "samplePointer", acqData.samplePointer)

    # fov and encoding size
    write(file, "fov", acqData.fov)
    write(file, "encodingSize", acqData.encodingSize)

    # kspace data
    rawdata_real = collect(reshape(reinterpret(Float64, vec(acqData.kdata)), (2,size(acqData.kdata)...)))
    write(file, "/rawdata", rawdata_real)

  end
end

function loadIBIFile(filename::AbstractString)
  acq = acquisitionData(MRIFileIBI(filename))
  return acq
end

# The follwing are helper functions that allow to store an entire Dict
# in an HDF5 file and load it.

export loadParams, saveParams

function saveParams(filename::AbstractString, path, params::Dict)
  h5open(filename, "w") do file
    saveParams(file, path, params)
  end
end

function saveParams(file, path, params::Dict)
  for (key,value) in params
    ppath = joinpath(path,string(key))
    if typeof(value) <: Bool
      write(file, ppath, UInt8(value))
      dset = file[ppath]
      attrs(dset)["isbool"] = "true"
    elseif typeof(value) <: Range
      write(file, ppath, [first(value),step(value),last(value)])
      dset = file[ppath]
      attrs(dset)["isrange"] = "true"
    elseif value == nothing
      write(file, ppath, "")
      dset = file[ppath]
      attrs(dset)["isnothing"] = "true"
    elseif typeof(value) <: Array{Any}
      write(file, ppath, [v for v in value])
    elseif typeof(value) <: Tuple
      write(file, ppath, [v for v in value])
      dset = file[ppath]
      attrs(dset)["istuple"] = "true"
    elseif typeof(value) <: AbstractLinearOperator
      continue
    else
      write(file, ppath, value)
    end
  end
end

function loadParams(filename::AbstractString, path)
  params = h5open(filename, "r") do file
   loadParams(file, path)
 end
  return params
end

function loadParams(file, path)
  params = Dict{Symbol,Any}()

  g = file[path]
  for obj in g
    key = last(splitdir(HDF5.name(obj)))
    data = read(obj)
    attr = attrs(obj)
    if exists(attr, "isbool")
      params[Symbol(key)] = Bool(data)
    elseif exists(attr, "isrange")
      if data[2] == 1
        params[Symbol(key)] = data[1]:data[3]
      else
        params[Symbol(key)] = data[1]:data[2]:data[3]
      end
    elseif exists(attr, "isnothing")
       params[Symbol(key)] = nothing
    elseif exists(attr, "istuple")
      params[Symbol(key)] = Tuple(data)
    else
      params[Symbol(key)] = data
    end
  end

  return params
end
