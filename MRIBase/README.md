# MRIBase

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRIBase/)

MRIBase defines the data types that the packages of MRIReco.jl exchange, such as `RawAcquisitionData`, `AcquisitionData` and `Trajectory`. It is one of the packages of [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl), does not depend on the others and is registered on its own:

```julia
using Pkg
Pkg.add("MRIBase")
```

The [documentation of MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/) covers all packages of the repository, and its [API reference for MRIBase](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRIBase/) lists the functions of this package. [Data Model](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/generated/explanations/data_model/) explains these data types.

Questions are welcome in the [discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions) and bug reports in the [issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues) of the repository. The [contributing guide](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/CONTRIBUTING.md) explains how to develop and test the packages.
