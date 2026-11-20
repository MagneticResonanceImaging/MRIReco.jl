# MRIBase

MRIBase defines the data types that the other packages build on: raw data organized like an ISMRMRD file (`RawAcquisitionData`), k-space data prepared for reconstruction (`AcquisitionData`) and k-space trajectories (`Trajectory`), together with functions to access and convert them. MRIBase does not depend on other packages of this repository. MRIReco re-exports it, so `using MRIReco` makes all names on this page available.

## Raw Data

```@docs
MRIBase.RawAcquisitionData
MRIBase.Profile
MRIBase.AcquisitionHeader
MRIBase.EncodingCounters
MRIBase.Limit
MRIBase.MeasurementDependency
MRIBase.CoilDescription
MRIBase.RawAcquisitionData(::AcquisitionData)
MRIBase.AcquisitionData(::RawAcquisitionData)
```

```@autodocs
Modules = [MRIBase]
Pages = ["Datatypes/RawAcqData.jl"]
Order = [:function]
Private = false
```

## Flags

The guide [Select and Prepare Profiles](@ref prepare-profiles) shows how to use the flags. Each name in
[`FLAGS`](@ref) is also exported as a `UInt64` bit mask, for example `ACQ_IS_NOISE_MEASUREMENT`, for working with
`head.flags` directly.

```@autodocs
Modules = [MRIBase]
Pages = ["Datatypes/Flags.jl"]
Private = false
```

## AcquisitionData

The guide [Create AcquisitionData from Arrays](@ref acqdata-from-arrays) shows how to create it from arrays, and [Data Model](@ref data-model) explains its fields.

```@autodocs
Modules = [MRIBase]
Pages = ["Datatypes/AcqData.jl"]
Private = false
```

## Trajectories

The guide [Create a Trajectory](@ref create-trajectory) shows how to create trajectories.

```@autodocs
Modules = [MRIBase]
Pages = ["Trajectories/Trajectories.jl"]
Private = false
```

### Predefined Trajectories

```@autodocs
Modules = [MRIBase]
Pages = ["2D/Cartesian2D.jl", "2D/EPI.jl", "2D/Radial2D.jl", "2D/Spiral2D.jl", "2D/Spiral2DVariableDens.jl",
         "2D/OneLine2D.jl", "2D/Spiral2DDualDens.jl", "2D/Spiral2DPerturbed.jl", "3D/Cartesian3D.jl",
         "3D/Kooshball.jl", "3D/StackofStars.jl"]
Private = false
```
