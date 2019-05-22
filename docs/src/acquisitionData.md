# Acquisition Data

There are two different forms of acquisition data types in `MRIReco`:
* `RawAcquisitionData`
* `AcquisitionData`
While the former is used to hold the data in the form, how it will be written
out from the scanner, the later has already performed some data permutations
bringing the data into the shape how the reconstruction expects it.

## Raw

The `RawAcquisitionData` is a data type that closely resembles the ISMRMRD data
format. It looks like
```julia
mutable struct RawAcquisitionData
  params::Dict{String, Any}
  profiles::Vector{Profile}
end
```
with
```julia
mutable struct Profile
  head::AcquisitionHeader
  traj::Array{Float32,2}
  data::Array{Complex{Float32},2}
end
```
The `params` member of `RawAcquisitionData` is basically the a flattened dictionary
derived from the XML part of an ISMRMRD file. A `Profile` describes the data
measured after a single excitation during an MRI experiment. It has members
`head`, `traj`, and `data`, which exactly resemble the structures specified
by the ISMRMRD file format. See [here](???) for more details.

## Preprocessed

The `RawAcquisitionData` can be preprocessed into a form, which makes it more
convenient for reconstruction algorithms. The `AcquisitionData` type looks like
```julia
mutable struct AcquisitionData
  sequenceInfo::Dict{Symbol,Any}
  traj::Vector{Trajectory}
  kdata::Array{Matrix{ComplexF64},3}
  subsampleIndices::Vector{Array{Int64}}
  encodingSize::Vector{Int64}
  fov::Vector{Float64}
end
```
It consists of the sequence informations stored in a dictionary, the k-space
trajectory, the k-space data, and several parameters describing the dimension of the data and some additional index vectors.

The k-space data `kdata` has three dimensions encoding
1. dim : contrasts/echoes
2. dim : slices
3. dim : repetitions
Each element is a matrix encoding
1. dim : k-space nodes
2. dim : channels/coils

In case of undersampled data, the subsampling indices are stored in `subsampleIndices`.
One check if the data is undersampled by checking if `isempty(subsampleIndices)`.

The encoded space is stored in the field `encodingSize`. It is especially relevant
for non-Cartesian trajectories where it is not clear upfront, how large the grid
size for reconstruction should be chosen. Finally `fov` describes the physical lengths
of the encoding grid.
