# MRIReco.jl

*Magnetic Resonance Imaging Reconstruction*

## Introduction

MPIReco is a Julia packet for magnetic resonance imaging. It contains algorithms for the simulation and reconstruction of MRT data and is both easy to use and flexibly expandable.

Both direct and iterative methods are available for image reconstruction. In particular, modern compressed sensing algorithms such as ADMM can be used.

The MRT imaging operator can be set up for a variety of scanning patterns (cartesian, spiral, radial, ...) and can take into account field inhomogeneity as well as the use of coil arrays. The operator can be quickly evaluated using NFFT-based methods.

One strength of the package is that it is strongly modular and uses high quality Julia packages. These are e.g.
 * NFFT.jl and FFTW.jl for fast Fourier transformations
 * Wavelets.jl for sparsification
 * LinearOperators.jl in order to be able to divide the imaging operator modularly into individual parts
 * RegularizedLeastSquares.jl for modern algorithms for solving linear optimization problems

This interaction allows new algorithms to be easily integrated into the software framework. It is not necessary to program in C/C++ but the advantages of the scientific high-level language Julia can be used.

!!! note
    MRIReco.jl is work in progress and in some parts not entirely optimized. In particular the FFT and NFFT implementation are currently limited to the CPU and do not support
    GPU acceleration yet.

## Installation

Start julia and open the package mode by entering `]`. Then enter
```julia
add MRIReco
```
This will install `MRIReco` and all its dependencies. If you want to develop
`MRIReco` itself you can checkout `MRIReco` by calling
```julia
add MRIReco
```
More information on how to develop a package can be found in the Julia documentation.

## Plotting

On purpose `MRIReco` is not depending on a particular plotting package since there
are various plotting packages in the Julia ecosystem. Within the examples outlined
in the tutorial we will use `PyPlot` for plotting but you may prefer using
the `Plots` package. You can add both packages the same way
`MRIReco` has been added.

## Tutorial

There is a Jupyter-based tutorial on `MRIReco` at

[https://github.com/MagneticResonanceImaging/MRIRecoTutorial](https://github.com/MagneticResonanceImaging/MRIRecoTutorial)

that has been presented at the ISMRM conference in Montreal 2019. Since the API
has slightly changed, we, however recommend that you read this documentation and
in particular execute the example scripts as is described in the [Getting Started](@ref) section.
