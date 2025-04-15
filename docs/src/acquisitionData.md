# Acquisition Data

There are two different forms of acquisition data types in `MRIReco`:
* `RawAcquisitionData`
* `AcquisitionData`
While the former is used to hold the data in the form, how it will be written
out from the scanner, the later has already performed some data permutations
bringing the data into the shape how the reconstruction expects it.

## Raw Data

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
The `params` member of `RawAcquisitionData` is basically a flattened dictionary
derived from the XML part of an ISMRMRD file. A `Profile` describes the data
measured after a single excitation during an MRI experiment. It has members
`head`, `traj`, and `data`, which exactly resemble the structures specified
by the ISMRMRD file format.

AcquisitionHeader has exactly the same structure as ISMRMRD. You can find more information about it [here](https://ismrmrd.readthedocs.io/en/latest/mrd_raw_data.html#acquisitionheader)

Two fields are especially important in it :
- flags
- idx

### flags

The flags field in the AcquisitionHeader is a 64 bit mask that can be used to indicate specific attributes of the corresponding readout. One usage of these flags is to reverse the signal during conversion from RawAcquisitionData to AcquisitionData if the flag "ACQ_IS_REVERSE" is set.

```julia
FLAGS = Dict(
    "ACQ_FIRST_IN_ENCODE_STEP1"                => 1,
    "ACQ_LAST_IN_ENCODE_STEP1"                 => 2,
    "ACQ_FIRST_IN_ENCODE_STEP2"                => 3,
    "ACQ_LAST_IN_ENCODE_STEP2"                 => 4,
    "ACQ_FIRST_IN_AVERAGE"                     => 5,
    "ACQ_LAST_IN_AVERAGE"                      => 6,
    "ACQ_FIRST_IN_SLICE"                       => 7,
    "ACQ_LAST_IN_SLICE"                        => 8,
    "ACQ_FIRST_IN_CONTRAST"                    => 9,
    "ACQ_LAST_IN_CONTRAST"                     => 10,
    "ACQ_FIRST_IN_PHASE"                       => 11,
    "ACQ_LAST_IN_PHASE"                        => 12,
    "ACQ_FIRST_IN_REPETITION"                  => 13,
    "ACQ_LAST_IN_REPETITION"                   => 14,
    "ACQ_FIRST_IN_SET"                         => 15,
    "ACQ_LAST_IN_SET"                          => 16,
    "ACQ_FIRST_IN_SEGMENT"                     => 17,
    "ACQ_LAST_IN_SEGMENT"                      => 18,
    "ACQ_IS_NOISE_MEASUREMENT"                 => 19,
    "ACQ_IS_PARALLEL_CALIBRATION"              => 20,
    "ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING"  => 21,
    "ACQ_IS_REVERSE"                           => 22,
    "ACQ_IS_NAVIGATION_DATA"                   => 23,
    "ACQ_IS_PHASECORR_DATA"                    => 24,
    "ACQ_LAST_IN_MEASUREMENT"                  => 25,
    "ACQ_IS_HPFEEDBACK_DATA"                   => 26,
    "ACQ_IS_DUMMYSCAN_DATA"                    => 27,
    "ACQ_IS_RTFEEDBACK_DATA"                   => 28,
    "ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA"    => 29,
    "ACQ_COMPRESSION1"                         => 53,
    "ACQ_COMPRESSION2"                         => 54,
    "ACQ_COMPRESSION3"                         => 55,
    "ACQ_COMPRESSION4"                         => 56,
    "ACQ_USER1"                                => 57,
    "ACQ_USER2"                                => 58,
    "ACQ_USER3"                                => 59,
    "ACQ_USER4"                                => 60,
    "ACQ_USER5"                                => 61,
    "ACQ_USER6"                                => 62,
    "ACQ_USER7"                                => 63,
    "ACQ_USER8"                                => 64
)
```

You can check the flags of a profile with `flags_of(p:Profile OR head::AcquisitionHeader)` or `flag_is_set` and manipulate them with thus functions :
- `flag_set!(obj::Profile, flag)`
- `flag_remove!(obj::Profile, flag)`
- `flag_remove_all!(obj::Profile)`


Thus functions can also been directly applied to the `AcquisitionHeader` and you can set! or remove! multiple flags at once :

```julia
flag_set!(head::AcquisitionHeader, [FLAGS["ACQ_USER8"], FLAGS["ACQ_IS_REVERSE"]])
# OR
flag_set!(head::AcquisitionHeader, [FLAGS["ACQ_USER8"],FLAGS["ACQ_IS_REVERSE"]])
```

Alternatively, you can set the flags directly with specific Const variable :
```julia
h = AcquisitionHeader()
h.flags = ACQ_IS_PARALLEL_CALIBRATION | ACQ_IS_NOISE_MEASUREMENT
flags_of(h)
h.flags = h.flags & ~ACQ_IS_NOISE_MEASUREMENT # remove the flag ACQ_IS_NOISE_MEASUREMENT
flags_of(h)
```

### idx

MR acquisitions often loop through a set of counters (e.g. phase encodes) in a complete experiment. The following encoding counters are referred to by the idx field in the AcquisitionHeader ([See the ISMRMRD documentation](https://ismrmrd.readthedocs.io/en/latest/mrd_raw_data.html#mrd-encodingcounters))

## Preprocessed Data

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

## Manipulating rawData

### Extract a subset of profiles with flags

A subset of rawData can be extracted according to the flags using the function `filter_raw_by_flags`. For example this function can be used to extract :
- Noise : `raw_noise = filter_raw_by_flags(raw,"ACQ_IS_NOISE_MEASUREMENT")`
- Extract calibration data : `raw_noise = filter_raw_by_flags(raw,["ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING","ACQ_IS_PARALLEL_CALIBRATION"])`

If you want to retrieve the data without the noise and calibration data you can use the following command :

`raw = remove_raw_by_flags(raw,["ACQ_IS_NOISE_MEASUREMENT","ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING","ACQ_IS_PARALLEL_CALIBRATION"])`

by default the following flags are removed :
```
"ACQ_IS_NOISE_MEASUREMENT", 
"ACQ_IS_PARALLEL_CALIBRATION", 
"ACQ_IS_NAVIGATION_DATA", 
"ACQ_IS_PHASECORR_DATA", 
"ACQ_IS_DUMMYSCAN_DATA", 
"ACQ_IS_PHASE_STABILIZATION_REFERENCE", "ACQ_IS_PHASE_STABILIZATION"
```

## Pre-processing tools

### Remove readout oversampling

On clinical MR scanners, the raw data are often oversampled in the readout direction. This oversampling is not necessary for image reconstruction and can be removed to save memory and computation time. This step is generally applied right after the acquisition of each line.
The current implementation is performed on the `RawAcquisitionData` (before conversion to `AcquisitionData`) and is based on the following parameters:
- `rawData.params["encodedFOV"][1]`
- `rawData.params["reconFOV"][1]`

If an oversampling factor of 2 is applied during acquisition, the `encodedFOV` will be 2 times larger than `reconFOV` along the readout dimension.

```
!!! note
    Be careful to remove all the profiles from a `RawAcquisitionData` that are not supposed to be of the same size; like noise acquisition or navigators. You can use the function `raw = remove_raw_by_flags(raw)`.
```