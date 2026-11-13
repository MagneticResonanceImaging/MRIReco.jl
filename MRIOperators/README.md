# MRIOperators

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRIOperators/)

MRIOperators implements the encoding operators of MRIReco.jl, which map an image to k-space data. Besides Fourier operators for Cartesian and non-Cartesian trajectories, it provides operators for coil sensitivities, off-resonance, relaxation and subspace models. It is one of the packages of [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl), builds on MRIBase and is registered on its own:

```julia
using Pkg
Pkg.add("MRIOperators")
```

The [documentation of MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/) covers all packages of the repository, and its [API reference for MRIOperators](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRIOperators/) lists the functions of this package. [Signal Model and Encoding Operators](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/explanations/signal_model/) explains these operators, and [GPU Acceleration](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/explanations/gpu/) describes how they run on GPUs.

Questions are welcome in the [discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions) and bug reports in the [issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues) of the repository. The [contributing guide](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/CONTRIBUTING.md) explains how to develop and test the packages.
