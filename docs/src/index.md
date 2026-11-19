```@raw html
---
layout: home

hero:
  name: MRIReco.jl
  tagline: Magnetic Resonance Imaging Reconstruction in Julia
  image:
    src: /logo.png
    alt: MRIReco.jl
  actions:
    - theme: brand
      text: Getting Started
      link: /gettingStarted
    - theme: alt
      text: Tutorials
      link: /generated/tutorials/direct
    - theme: alt
      text: View on GitHub
      link: https://github.com/MagneticResonanceImaging/MRIReco.jl

features:
  - title: Raw Data and Files
    details: Read and write ISMRMRD files, read Bruker datasets and prepare raw data for reconstruction.
    link: /generated/howto/load_save
  - title: Iterative Reconstruction
    details: Reconstruct undersampled data with the solvers and regularization terms of RegularizedLeastSquares.jl.
    link: /generated/tutorials/iterative
  - title: Parallel Imaging
    details: Estimate coil sensitivities with ESPIRiT and reconstruct one image from the data of all coils, also combined with compressed sensing.
    link: /generated/tutorials/multicoil
  - title: GPU Acceleration
    details: Run reconstructions on a GPU by setting the reconstruction parameter <code>arrayType</code>.
    link: /generated/howto/gpu
---
```

MRIReco.jl is a framework for MRI image reconstruction written in Julia. It provides data types for raw data and for k-space data prepared for reconstruction, encoding operators, and direct as well as iterative reconstruction methods, including parallel imaging and compressed sensing. Further packages read and write raw data files, generate sampling patterns, estimate coil sensitivities and simulate data.

## Packages

MRIReco.jl consists of several packages that are developed in the same repository and registered individually. The last column lists the packages of this repository that each package depends on. MRIReco re-exports its dependencies MRIBase and MRIOperators, so their names are available after `using MRIReco`; the other packages are independent of MRIReco and are installed and loaded separately.

| Package | Purpose | Depends on |
|:---|:---|:---|
| [MRIReco](API/MRIReco.md) | Image reconstruction | MRIBase, MRIOperators |
| [MRIBase](API/MRIBase.md) | Data types for raw data, acquisition data and trajectories | — |
| [MRIOperators](API/MRIOperators.md) | Encoding operators | MRIBase |
| [MRIFiles](API/MRIFiles.md) | Reading and writing raw data files | MRIBase |
| [MRISampling](API/MRISampling.md) | Undersampling patterns | MRIBase |
| [MRICoilSensitivities](API/MRICoilSensitivities.md) | Coil sensitivity estimation and coil compression | MRIBase |
| [MRISimulation](API/MRISimulation.md) | Simulation of MRI data | MRIBase, MRIOperators |

## How the Documentation Is Organized

- [Getting Started](@ref) shows how to install the packages and run a first reconstruction.
- The tutorials walk through complete reconstructions step by step.
- The how-to guides show how to solve specific tasks, grouped by package.
- The explanations describe the design of the packages and the concepts behind them.
- The reference lists the API of each package.

## Citing MRIReco.jl

If you use MRIReco.jl in your work, please cite

> T. Knopp and M. Grosser, "MRIReco.jl: An MRI reconstruction framework written in Julia", Magnetic Resonance in Medicine 86, 1633–1646 (2021), [doi:10.1002/mrm.28792](https://doi.org/10.1002/mrm.28792).

A BibTeX entry is in [`MRIReco.bib`](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/MRIReco.bib).

## Getting Help

Questions are welcome in [GitHub Discussions](https://github.com/MagneticResonanceImaging/MRIReco.jl/discussions), and bugs can be reported as [issues](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues). To contribute changes, see the [contributing guide](contributing.md).
