# MRICoilSensitivities

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRICoilSensitivities/)

MRICoilSensitivities estimates coil sensitivities, for example with ESPIRiT, and compresses the data of many coils into fewer virtual coils. It is one of the packages of [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl), builds on MRIBase and is registered on its own:

```julia
using Pkg
Pkg.add("MRICoilSensitivities")
```

The [documentation of MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/) covers all packages of the repository, and its [API reference for MRICoilSensitivities](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRICoilSensitivities/) lists the functions of this package. [Multi-Coil Reconstruction](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/generated/tutorials/multicoil/) uses it to estimate the sensitivities.

Questions are welcome in the [discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions) and bug reports in the [issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues) of the repository. The [contributing guide](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/CONTRIBUTING.md) explains how to develop and test the packages.
