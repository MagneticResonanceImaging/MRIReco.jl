export CFLFile, data

struct CFLFile
  filenameBase::String
end

#
# Read in recon data stored in filenameBase.cfl (complex float)
# based on dimensions stored in filenameBase.hdr.
#
function data(f::CFLFile)

    dims = readCFLHeader(f.filenameBase)

    fid = open(f.filenameBase*".cfl")

    dataFloat = read!(fid, Array{Float32}(undef,2*prod(dims)))
    data = zeros(ComplexF64, dims...)
    data[:] = dataFloat[1:2:end] + 1im*dataFloat[2:2:end]

    close(fid)
    return data
end

#
# reader header file and return size of the data
#
function readCFLHeader(filenameBase::String)
    fid = open(filenameBase*".hdr","r")

    line = getNextLine(fid)
    dims = map(x->parse(Int32,x), split(line))

    close(fid)
    return dims
end

#
# return the next line in <fid>, which is not a comment
#
function getNextLine(fid)
  line = lstrip(readuntil(fid, "\n"))
  while(line[1] == '#')
    line = readuntil(fid, "\n")
  end

  return line
end
