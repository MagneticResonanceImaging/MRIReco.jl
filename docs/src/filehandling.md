# File Handling

TODO

MRI [Acquisition Data](@ref) can not only be generated from [Simulation](@ref)
but also from files. Currently, MRIReco supports the ISMRMRD data format and
a custom MRILib specific data format. Both formats are a subtype of `MRIFile` and
implement the functions
* `trajectory(f::MRIFile)`
* `sequence(f::MRIFile)`
* `rawdata(f::MRIFile)`
* `AcquisitionData(f::MRIFile)`
which allow to access specific parts of the MRIFile. The last function `AcquisitionData`
returns an `AcquisitionData` data object which can be used directly in the
[Reconstruction](@ref) methods.
