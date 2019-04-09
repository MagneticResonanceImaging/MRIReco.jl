export trajectory

export DFFile

struct DFFile <: MRIFile
  datafilename::String
  trajfilename::String
end

# reading

function trajectory(f::DFFile; kargs...)
  numSamplingPerProfile, numProfiles, nodes = open(f.trajfilename,"r") do fd
    # tmp1,numSamplingPerProfile,numProfiles,tmp2= read(fd,Int32,4)
    tmp1,numSamplingPerProfile,numProfiles,tmp2= read!(fd,Array{Int32}(undef,4))
    # nodes = read(fd,Float32,2,numSamplingPerProfile,numProfiles)
    nodes = read!(fd, Array{Float32}(undef, 2, numSamplingPerProfile,numProfiles))
    return numSamplingPerProfile, numProfiles, nodes
  end

  return Trajectory(reshape(nodes,2,:), numProfiles, numSamplingPerProfile; kargs...)
end

function rawdata(f::DFFile)
  data = open(f.datafilename,"r") do fd
    # tmp1,numSamplingPerProfile,numProfiles,numSlices= read(fd,Int32,4)
    tmp1,numSamplingPerProfile,numProfiles,numSlices= read!(fd,Array{Int32}(undef,4))
    # data = read(fd,ComplexF32,numSamplingPerProfile,numProfiles,numSlices)
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
