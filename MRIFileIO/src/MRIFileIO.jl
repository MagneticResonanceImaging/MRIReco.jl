module MRIFileIO

using Graphics: @mustimplement
using Reexport

using MRIBase
using Reexport
using FileIO
using LinearAlgebra
using HDF5
using LightXML
using NIfTI


# With the current design, the following type is not necessary, to be discussed
abstract type MRIFile end

include("Parameters.jl")
include("ISMRMRD.jl")
include("DFFile.jl")
include("CFLFile.jl")
#include("Nifti.jl")
include("Brukerfile.jl")

end # module