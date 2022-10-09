module MRIFiles

using Reexport

using MRIBase
using MRIBase.NFFTTools.AbstractNFFTs: @mustimplement # magic line
using Reexport
using FileIO
using LinearAlgebra
using HDF5
using LightXML

abstract type MRIFile end
@mustimplement MRIBase.RawAcquisitionData(::MRIFile)

include("ISMRMRD/ISMRMRD.jl")
include("Bruker/Bruker.jl")
include("Siemens/Siemens.jl")
include("Philips/Philips.jl")
include("GE/GE.jl")
include("Others/Others.jl")

#include("Nifti.jl")

end # module