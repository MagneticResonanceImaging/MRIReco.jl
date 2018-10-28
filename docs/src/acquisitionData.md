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
  idx::Array{Int64}
end
```
The composite type thus consists of the imaging sequence, the k-space data,
several parameters describing the dimension of the data and some additional
index vectors.
