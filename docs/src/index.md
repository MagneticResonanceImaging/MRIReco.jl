# MRIReco.jl

*Magnetic Resonance Imaging Reconstruction in Julia*

MRIReco.jl is a framework for MRI image reconstruction written in Julia. It provides data types for raw data and for k-space data prepared for reconstruction, encoding operators, and direct as well as iterative reconstruction methods, including parallel imaging and compressed sensing. Further packages read and write raw data files, generate sampling patterns, estimate coil sensitivities and simulate data.

## Packages

MRIReco.jl consists of several packages that are developed in the same repository and registered individually. `using MRIReco` also loads MRIBase and MRIOperators; the other packages are installed and loaded separately.

| Package | Purpose |
|:---|:---|
| [MRIReco](API/MRIReco.md) | Image reconstruction |
| [MRIBase](API/MRIBase.md) | Data types for raw data, acquisition data and trajectories |
| [MRIOperators](API/MRIOperators.md) | Encoding operators |
| [MRIFiles](API/MRIFiles.md) | Reading and writing raw data files |
| [MRISampling](API/MRISampling.md) | Undersampling patterns |
| [MRICoilSensitivities](API/MRICoilSensitivities.md) | Coil sensitivity estimation and coil compression |
| [MRISimulation](API/MRISimulation.md) | Simulation of MRI data |

## How the Documentation Is Organized

- [Getting Started](@ref) shows how to install the packages and run a first reconstruction.
- The tutorials walk through complete reconstructions step by step.
- The how-to guides show how to solve specific tasks, grouped by package.
- The reference lists the API of each package.

The pages under "Legacy Pages" come from the previous version of this documentation and are being rewritten.

## Citing MRIReco.jl

If you use MRIReco.jl in your work, please cite

> T. Knopp and M. Grosser, "MRIReco.jl: An MRI reconstruction framework written in Julia", Magnetic Resonance in Medicine 86, 1633–1646 (2021), [doi:10.1002/mrm.28792](https://doi.org/10.1002/mrm.28792).

A BibTeX entry is in [`MRIReco.bib`](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/MRIReco.bib).

## Getting Help

Questions are welcome in [GitHub Discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions), and bugs can be reported as [issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues).
