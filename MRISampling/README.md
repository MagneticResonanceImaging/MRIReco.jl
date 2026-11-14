# MRISampling

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRISampling/)

MRISampling generates undersampling patterns and applies them to acquisition data. It is one of the packages of [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl), builds on MRIBase and is registered on its own:

```julia
using Pkg
Pkg.add("MRISampling")
```

The [documentation of MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/) covers all packages of the repository, and its [API reference for MRISampling](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRISampling/) lists the functions of this package. [Iterative Reconstruction](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/generated/tutorials/iterative/) uses it to undersample data.

Questions are welcome in the [discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions) and bug reports in the [issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues) of the repository. The [contributing guide](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/CONTRIBUTING.md) explains how to develop and test the packages.
