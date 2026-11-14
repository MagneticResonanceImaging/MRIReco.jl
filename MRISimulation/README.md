# MRISimulation

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRISimulation/)

MRISimulation simulates k-space data of an image, taking into account coil sensitivities, off-resonance and relaxation, with the Fourier encoding operators of MRIOperators. It is one of the packages of [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl), builds on MRIBase and MRIOperators and is registered on its own:

```julia
using Pkg
Pkg.add("MRISimulation")
```

The [documentation of MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/) covers all packages of the repository, and its [API reference for MRISimulation](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRISimulation/) lists the functions of this package.

Questions are welcome in the [discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions) and bug reports in the [issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues) of the repository. The [contributing guide](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/CONTRIBUTING.md) explains how to develop and test the packages.
