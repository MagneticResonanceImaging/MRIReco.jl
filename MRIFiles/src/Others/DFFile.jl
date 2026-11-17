
export DFFile

"""
    DFFile(datafilename::String, trajfilename::String)

reference to 2D non-Cartesian data stored in two binary files, `trajfilename` with the k-space nodes and
`datafilename` with the samples. Both files start with four `Int32` values. In the trajectory file, the second and
third are the number of samples per profile and the number of profiles, followed by the nodes as `Float32` values
(2 × samples × profiles). In the data file, the second to fourth are the numbers of samples per profile, profiles and
slices, followed by the samples as `ComplexF32` values.

`trajectory(f)` returns the trajectory, and keyword arguments are passed on to `Trajectory`. `MRIFiles.rawdata(f)`
returns the samples as an `Array{ComplexF64,3}`, and `MRIFiles.rawdata(f, slice)` those of one slice.
"""
struct DFFile <: MRIFile
  datafilename::String
  trajfilename::String
end

# reading

function MRIBase.trajectory(f::DFFile; kargs...)
  numSamplingPerProfile, numProfiles, nodes = open(f.trajfilename,"r") do fd
    tmp1, numSamplingPerProfile, numProfiles, tmp2 = Int64.(read!(fd,Array{Int32}(undef,4)))
    nodes = read!(fd, Array{Float32}(undef, 2, numSamplingPerProfile,numProfiles))
    return numSamplingPerProfile, numProfiles, nodes
  end

  return Trajectory(reshape(nodes,2,:), numProfiles, numSamplingPerProfile; kargs...)
end

function rawdata(f::DFFile)
  data = open(f.datafilename,"r") do fd
    tmp1, numSamplingPerProfile, numProfiles, numSlices = Int64.(read!(fd,Array{Int32}(undef,4)))
    data = read!(fd, Array{ComplexF32}(undef,numSamplingPerProfile,numProfiles,numSlices))
    return data
  end

  return convert(Array{ComplexF64},data)
end

function rawdata(f::DFFile, slice::Int)
  data = open(f.datafilename,"r") do fd
    # tmp1,numSamplingPerProfile,numProfiles,numSlices= read(fd,Int32,4)
    tmp1,numSamplingPerProfile,numProfiles,numSlices= read!(fd,Array{Int32}(undef,4))
    # data = read(fd,ComplexF32,numSamplingPerProfile,numProfiles,numSlices)
    data = read!(fd, Array{ComplexF32}(undef,numSamplingPerProfile,numProfiles,numSlices))
    return data
  end

  return convert(Array{ComplexF64}, data[:,:,slice])
end
