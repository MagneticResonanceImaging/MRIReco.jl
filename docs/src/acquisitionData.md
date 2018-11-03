# Acquisition Data

All acquisition data is stored in the a type that looks like this
```julia
mutable struct AcquisitionData{S<:AbstractSequence}
  seq::S
  kdata::Vector{ComplexF64}
  numEchoes::Int64
  numCoils::Int64
  numSlices::Int64
  samplePointer::Vector{Int64}
  subsampleIndices::Array{Int64}
  encodingSize::Vector{Int64}
  fov::Vector{Float64}
end
```
The composite type consists of the imaging sequence, the k-space data,
several parameters describing the dimension of the data and some additional
index vectors.

The k-space data `kdata` is flattened into a 1D vector but it represents data
from a 4D space with dimensions
1. kspace nodes
2. echo times
3. coils
4. slices / repetitions
The reason to use a flattened 1D data is that the number k-space nodes needs not
to be constant for different echo times. The entry point to the data is stored
in the index vector `samplePointer`. It has length
```julia
numEchoes * numCoils * numSlices
```
and gives for each combination of echo, coil and slice the corresponding index,
where the k-space data starts. The end-point can be obtained by incrementing the index
by one.

In case of undersampled data, the subsampling indices are stored in `subsampleIndices`.
One check if the data is undersampled by checking if `isempty(subsampleIndices)`.

The encoded space is stored in the field `encodingSize`. It is especially relevant
for non-Cartesian trajectories where it is not clear upfront, how large the grid
size for reconstruction can be chosen. Finally `fov` describes the physical lengths
of the encoding grid.
