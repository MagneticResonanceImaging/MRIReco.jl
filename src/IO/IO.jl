# The purpose of the abstract MRI file is to
# implement different MRI file formats, in particular the ISMRMRD format
# Since the later is a little bit more complex I will start with a homegrown
# stuff that we can replace later.

export MRIFile

abstract type MRIFile end


include("ISMRMRD.jl")
include("MRIFileIBI.jl")
include("RecoFileIBI.jl")
include("DFFile.jl")
include("CFLFile.jl")
include("Nifti.jl")
include("ImageExport.jl")
