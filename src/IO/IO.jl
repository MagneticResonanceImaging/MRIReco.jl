export MRIFile

# With the current design, the following type is not necessary, to be discussed
abstract type MRIFile end

include("ISMRMRD.jl")
include("DFFile.jl")
include("CFLFile.jl")
include("Nifti.jl")
include("Brukerfile.jl")
