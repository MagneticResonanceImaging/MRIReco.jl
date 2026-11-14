# MRIReco

![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build status](https://github.com/MagneticResonanceImaging/MRIReco.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/MagneticResonanceImaging/MRIReco.jl/actions/workflows/ci.yml)
[![codecov.io](http://codecov.io/github/MagneticResonanceImaging/MRIReco.jl/coverage.svg?branch=master)](http://codecov.io/github/MagneticResonanceImaging/MRIReco.jl?branch=master)


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://magneticresonanceimaging.github.io/MRIReco.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://magneticresonanceimaging.github.io/MRIReco.jl/dev)

MRIReco.jl is a framework for MRI image reconstruction written in Julia. It provides data types for raw data and for k-space data prepared for reconstruction, encoding operators, and direct as well as iterative reconstruction methods, including parallel imaging and compressed sensing.

## Packages

The repository contains several packages, which are registered individually. The last column lists the packages of this repository that each package depends on.

| Package | Purpose | Depends on |
|:---|:---|:---|
| [MRIReco](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/API/MRIReco/) | Image reconstruction | MRIBase, MRIOperators |
| [MRIBase](https://github.com/MagneticResonanceImaging/MRIReco.jl/tree/master/MRIBase) | Data types for raw data, acquisition data and trajectories | — |
| [MRIOperators](https://github.com/MagneticResonanceImaging/MRIReco.jl/tree/master/MRIOperators) | Encoding operators | MRIBase |
| [MRIFiles](https://github.com/MagneticResonanceImaging/MRIReco.jl/tree/master/MRIFiles) | Reading and writing raw data files | MRIBase |
| [MRISampling](https://github.com/MagneticResonanceImaging/MRIReco.jl/tree/master/MRISampling) | Undersampling patterns | MRIBase |
| [MRICoilSensitivities](https://github.com/MagneticResonanceImaging/MRIReco.jl/tree/master/MRICoilSensitivities) | Coil sensitivity estimation and coil compression | MRIBase |
| [MRISimulation](https://github.com/MagneticResonanceImaging/MRIReco.jl/tree/master/MRISimulation) | Simulation of MRI data | MRIBase, MRIOperators |

## Installation

```julia
using Pkg
Pkg.add("MRIReco")
```

`using MRIReco` also loads MRIBase and MRIOperators. The other packages are added and loaded separately, for example MRIFiles to read raw data files. [Getting Started](https://magneticresonanceimaging.github.io/MRIReco.jl/dev/gettingStarted/) runs a first reconstruction.

## How to give credit
If you use this package please acknowledge us by citing:

T. Knopp and M. Grosser (2021). MRIReco.jl: An MRI Reconstruction Framework written in Julia. Magnetic Resonance in Medicine. 2021 [pdf](https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.28792) [url](https://doi.org/10.1002/mrm.28792) [*arXiv:2101.12624*](https://arxiv.org/abs/2101.12624)

A BibTeX file is directly contained within this package `MRIReco.bib`.

# Community Standards

This project is part of the Julia community and follows the [Julia community standards](https://julialang.org/community/standards/). 

The [contributing guide](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/CONTRIBUTING.md) explains how to report problems and how to develop and test the packages, and the [code of conduct](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/CODE_OF_CONDUCT.md) applies to everyone who takes part.
