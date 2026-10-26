# MRIBase

MRIBase defines the data types that the other packages build on: raw data organized like an ISMRMRD file (`RawAcquisitionData`), k-space data prepared for reconstruction (`AcquisitionData`) and k-space trajectories (`Trajectory`), together with functions to access and convert them. MRIReco re-exports MRIBase, so `using MRIReco` makes all names on this page available.

## Raw Data

```@docs
MRIBase.RawAcquisitionData
MRIBase.Profile
MRIBase.AcquisitionHeader
MRIBase.EncodingCounters
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

```@autodocs
Modules = [MRIBase]
Pages = ["Datatypes/Flags.jl"]
Private = false
```

## AcquisitionData

```@autodocs
Modules = [MRIBase]
Pages = ["Datatypes/AcqData.jl"]
Private = false
```

## Trajectories

```@autodocs
Modules = [MRIBase]
Pages = ["Trajectories/Trajectories.jl"]
Private = false
```

### Predefined Trajectories

```@autodocs
Modules = [MRIBase]
Pages = ["2D/Cartesian2D.jl", "2D/EPI.jl", "2D/Radial2D.jl", "2D/Spiral2D.jl", "2D/Spiral2DVariableDens.jl",
         "3D/Cartesian3D.jl", "3D/Kooshball.jl", "3D/StackofStars.jl"]
Private = false
```
